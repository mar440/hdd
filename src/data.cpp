#include "../include/data.hpp"
#include "../include/element.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/stiffnessMatrix.hpp"

#include <vtkCell.h>
#include <string>



using namespace std;
using namespace Eigen;

namespace pt = boost::property_tree;


Data::Data(MPI_Comm* _pcomm): m_domain(_pcomm)
{
  m_pcomm = _pcomm;
  m_p_interfaceOperatorB = nullptr;

  MPI_Comm_rank(*m_pcomm, &m_mpiRank);
  MPI_Comm_size(*m_pcomm, &m_mpiSize);

  m_container.resize(0);
  m_defectPerSubdomains.resize(0);
  m_DirichletGlbDofs.resize(0); 


  std::string fname = "out_" + std::to_string(m_mpiRank) + ".txt";


  m_filestr.open (fname);
  m_p_backup = std::cout.rdbuf(); // back up cout's streambuf

  m_p_sbuf = m_filestr.rdbuf();   // get file's streambuf
  std::cout.rdbuf(m_p_sbuf);
  std::cout << "++++++++++++++++++++++++\n";
  std::cout << "SUBDOMAIN id. " << m_mpiRank << "\n";
  std::cout << "++++++++++++++++++++++++\n\n";
}

Data::~Data()
{
  if (m_p_interfaceOperatorB) delete m_p_interfaceOperatorB;
}

void Data::SymbolicAssembling(std::vector<int>& elemntIds)
{
  // gather global DOFs of subdomain elemnts to create
  // mapping vector
  m_container.insert(m_container.end(), elemntIds.begin(),
      elemntIds.end());
}

void Data::Finalize()
{
  std::cout.rdbuf(m_p_backup);        // restore cout's original streambuf

  m_filestr.close();
}

void Data::FinalizeSymbolicAssembling()
{
  // maping DOF: local to global
  linalg::unique(m_container);
  m_domain.SetMappingLoc2Glb(m_container);
  m_container.clear();
  m_container.shrink_to_fit();
  //
  m_domain.SetInterfaces();
  m_p_interfaceOperatorB = new InterfaceOperatorB(&m_domain);


  if (m_verboseLevel>0) m_dbg_printStats();


  std::string precondType =
    m_root.get<std::string>("solver.preconditioner");

  PRECONDITIONER_TYPE _precondType;

  if (precondType == "Dirichlet")
    _precondType = DIRICHLET;
  else if (precondType == "Lumped")
    _precondType = LUMPED;
  else if (precondType == "none")
    _precondType = NONE;
  else
  {
    _precondType = NONE;
    std::runtime_error("Unknown preconditioner type");
  }

  m_domain.InitStiffnessMatrix(_precondType);

}



void Data::NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS)
{
  m_domain.GetStiffnessMatrix()->AddElementContribution(glbIds, valLocK);
  m_domain.GetStiffnessMatrix()->AddRHSContribution(glbIds, valLocRHS);
}


void Data::FinalizeNumericAssembling()
{

  // DIRICHLET BOUNDARY CONDTION
  if (m_DirichletGlbDofs.size() == 0)
    std::runtime_error("Dirichlet BC must be provided \
        before \"FinalizeNumericAssembling\" is called");

  m_domain.SetDirichletDOFs(m_DirichletGlbDofs);

  // ASSEMBLE & FACTORIZE LINEAR OPERATOR
  m_domain.GetStiffnessMatrix()->FinalizeNumericPart(m_domain.GetDirichletDOFs());

  // DIRICHLET PRECONDITIONER
  m_domain.HandlePreconditioning();


  // NUMEBRING FOR GtG MATRIX
  m_SetKernelNumbering();
  // FETI COARSE SPACE 
  auto mpB = m_p_interfaceOperatorB;
  mpB->FetiCoarseSpace(_GetDefectPerSubdomains());








}



void Data::SetDirichletDOFs(std::vector<int>&v)
{
  //
  int _size = (int)v.size();
  m_DirichletGlbDofs.resize(_size); 
  //
  for (int dof = 0; dof < _size; dof++)
    m_DirichletGlbDofs[dof] = v[dof];


}



void Data::m_SetKernelNumbering()
{
  m_defectPerSubdomains.resize(m_mpiSize,-1);

  m_container.resize(m_mpiSize,m_domain.GetStiffnessMatrix()->GetDefect());


if (m_verboseLevel>2) {
  for (auto& iw : m_container) std::cout<< iw << ' ';
  std::cout << '\n';
}
  MPI_Alltoall(
      m_container.data(),1,MPI_INT,
      m_defectPerSubdomains.data(),1,MPI_INT,
      *m_pcomm);

if (m_verboseLevel>2) {
  for (auto& iw : m_defectPerSubdomains) std::cout<< iw << ' ';
  std::cout << '\n';
}

  m_container.resize(0);
  m_container.shrink_to_fit();

}

void Data::Solve(Eigen::VectorXd& solution)
{

  bool solverState = m_solver.pcpg(*this,solution);


  if (solverState)
    std::cout << "successful\n";
  else
    std::cout << "Iterative solver didn't finish successfully.\n";

  int dumpMatrices = m_root.get<int>("outputs.dumpMatrices",0);

  if (dumpMatrices!=0) _dumpTxtFiles(solution);

}



void Data::_dumpTxtFiles(Eigen::VectorXd& solution)
{
  m_domain.GetStiffnessMatrix()->PrintStiffnessMatrix("K_");
  m_domain.GetStiffnessMatrix()->PrintKernel("R_");
  m_domain.GetStiffnessMatrix()->PrintRHS("f_");
  m_domain.GetStiffnessMatrix()->PrintNullPivots("nullPivots_");
  m_p_interfaceOperatorB->printInterfaceDOFs("B_");
  m_p_interfaceOperatorB->printNeighboursRanks("neighbRanks_");

  std::string fname = "solution_" + std::to_string(m_mpiRank) + ".txt";
  tools::printMatrix(solution,fname);
}





///////////////////////////////////////
// build-in assembler for fast testing
///////////////////////////////////////



//void Data::assembly_elasticity(Mesh *mesh)
//{
//
//  int nel = mesh->getGlobalMesh()->GetNumberOfCells();
//  cout << "nel = " << nel << endl;
//  int dim = 2;
//  //double mat_E_mu[3] = {2.1e5, 0.3, 0.78500e9};
//  double mat_E_mu[3] = {1.0, 0.3, 1.0};
//
//  cout << "K assembly elasticity ..." << endl;
//  int nDOFs(0);
//  nDOFs = mesh->getGlobalMesh()->GetNumberOfPoints() * dim;
//  m_rhs.resize(nDOFs);
//  m_rhs.setZero();
//  std::vector<T> trK;
//  MatrixXd K_loc;
//  VectorXd f_loc;
//
//  for (int iCell = 0; iCell < nel; iCell++)
//  {
//    howMuchAssembled(iCell,nel);
//
//    auto cell = mesh->getGlobalMesh()->GetCell(iCell);
//    int nP = cell->GetNumberOfPoints();
//
//    Element *element;
//
//    switch (cell->GetCellType()){
//      case VTK_QUADRATIC_QUAD:
//        element = new QUADRATIC_QUAD;
//        break;
//      default:
//        continue;
//    }
//
//
//    element->assembly_elasticity(K_loc, f_loc, cell, mat_E_mu);
//
//    //    auto singVals = svd0(K_loc);
//    //    for (int vls = 0; vls < singVals.size(); vls++)
//    //      cout << singVals(vls) << ' ';
//    //    cout << endl;
//    //TODO the numbering for tetra 10 is (probably) reversed ...
//
//    vtkIdList *ids = cell->GetPointIds();
//
//
//    int i_glb, j_glb; //, i_loc, //j_loc;
//    for (int iDim = 0; iDim < dim; iDim++){
//      for (int jDim = 0; jDim < dim; jDim++){
//        for (int i = 0; i < nP; i++){
//          i_glb = dim * ids->GetId(i) + iDim;
////          i_loc = i_glb;//mesh->g2l_volume[i_glb];
//          if (jDim == 0)
//            m_rhs(i_glb) += f_loc(i + iDim * nP);
//          for (int j = 0; j < nP; j++){
//            j_glb = dim * ids->GetId(j) + jDim;
//            //j_loc = j_glb;// mesh->g2l_volume[j_glb];
//            trK.push_back(T(i_glb, j_glb, K_loc(i + iDim * nP, j + jDim * nP)));
//          }
//        }
//      }
//    }
//    delete element;
//  }
//
//
//  // global stiffness matrix for elasticity - initialization
//  m_K.resize(nDOFs, nDOFs);
//  m_K.setFromTriplets(trK.begin(), trK.end());
//
//
//
//  cout << " elasticity matrix assembled ... " << endl;
//}
//
//void Data::setDirichlet(Eigen::VectorXd& v,
//    std::vector<int> dirInd)
//{
//
//  for (int iD = 0; iD < static_cast<int>(dirInd.size()); iD++)
//    v[dirInd[iD]] = 0;
//
//}
//
//
//void Data::setDirichlet(
//    SpMat &mat,
//    std::vector<int> dirInd)
//{
//
//  VectorXd diagMat = mat.diagonal();
//  double meanVal = 1;// diagMat.sum() / diagMat.size();
//
//  for (int iD = 0; iD < static_cast<int>(dirInd.size()); iD++){
//
//    int kD = dirInd[iD];
//    for (SpMat::InnerIterator it(mat, kD); it; ++it)
//    {
//      if (it.row() == it.col())
//        it.valueRef() = meanVal;
//      else
//        it.valueRef() = 0;
//    }
//  }
//}


void Data::howMuchAssembled(int iCell, int nel){

  double aaa = static_cast<double> (100 * (iCell + 1)) / nel;
  if ((fmod(aaa, 5)) < 1e-2)
    std::cout << aaa << "% " << "(" << iCell << "," << nel <<
      ") ... fmod(aaa,5)" << (fmod(aaa, 5)) << std::endl;

}

void Data::m_dbg_printStats()
{
// test
  auto mpIO = m_p_interfaceOperatorB;

  int neqDual = m_domain.GetNumberOfDualDOFs();
  int neqPrimal = m_domain.GetNumberOfPrimalDOFs();

  int nRHS = 3;

  Eigen::MatrixXd vecD(neqDual,nRHS);
  for (int col = 0; col < nRHS; col++)
    vecD.col(col) = 10.0 * (col + 1) * Eigen::VectorXd::Ones(neqDual);
  Eigen::MatrixXd vecP(neqPrimal,nRHS);

  std::cout << "dual\n";
  for (int i = 0 ; i < neqDual; i++)
  {
    for (int j = 0 ; j < nRHS; j++)
    {
      std::cout << vecD(i,j) << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << " *" << std::endl;



  mpIO->multBt(vecD,vecP);

  std::cout << "primal\n";
  for (int i = 0 ; i < neqPrimal; i++)
  {
    for (int j = 0 ; j < nRHS; j++)
    {
      std::cout << vecP(i) << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << " *" << std::endl;

  mpIO->multB(vecP,vecD);
  
  std::cout << "dual\n";
  for (int i = 0 ; i < neqDual; i++)
  {
    for (int j = 0 ; j < nRHS; j++)
    {
      std::cout << vecD(i,j) << ' ';
    }
    std::cout << std::endl;
  }
  std::cout << " *" << std::endl;


}

//{
//    "solver" :
//    {
//        "preconditioner" : "Dirichlet",
//        "stopingCriteria" : 1e-4,
//        "maxNumbIter" : 200
//    },
//    "meshgenerator" :
//    {
//        "numberOfElements_x" : 10,
//        "numberOfElements_y" : 10,
//        "numberOfDomains_x" : 3,
//        "numberOfDomains_y" : 3,
//        "lenght_x" : 3,
//        "lenght_y" : 3,
//        "numberOfLevels" : 3
//    }
//}
//  std::cout << "PJ: in parser" << std::endl;
//
//// Short alias for this namespace
//
//// Create a root
//  pt::ptree m_root;
//
//// Load the json file in this ptree
//  pt::read_json(path2file, root);
//
//
//  // Read values
//  int height = root.get<int>("height", 0);
//  std::cout << "PJ: height: " << height << std::endl;
//  // You can also go through nested nodes
//  std::string msg = root.get<std::string>("some.complex.path");
//  std::cout << "PJ: msg: " << msg<< std::endl;



void Data::ParseJsonFile(std::string path2file)
{

  std::cout << "parsing \"hddConf.json\" file\n";
  pt::read_json(path2file, m_root);

  m_verboseLevel = m_root.get<int>("outputs.verbose",0);

//  std::string precond = m_root.get<std::string>("solver.preconditioner");
//  std::cout << "preconditioner: " << precond << std::endl;
//  int nx = m_root.get<int>("meshgenerator.numberOfElements_x",0);
//  std::cout << "meshgenerator.numberOfElements_x: " << nx << std::endl;

}
