#include "../include/data.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/stiffnessMatrix.hpp"
#include <string>



using namespace std;
using namespace Eigen;

namespace pt = boost::property_tree;


Data::Data(MPI_Comm* _pcomm): m_domain(_pcomm)
{
  m_comm = *_pcomm;
  m_p_interfaceOperatorB = nullptr;

  MPI_Comm_rank(m_comm, &m_mpiRank);

  m_container.resize(0);
  m_defectPerSubdomains.resize(0);
  m_DirichletGlbDofs.resize(0);


  std::string fname = "out_" + std::to_string(m_mpiRank) + ".txt";


  m_filestr.open(fname);
  m_p_backup = std::cout.rdbuf(); // back up cout's streambuf

  m_p_sbuf = m_filestr.rdbuf();   // get file's streambuf
  std::cout.rdbuf(m_p_sbuf);
  std::cout << "++++++++++++++++++++++++\n";
  std::cout << "SUBDOMAIN id. " << m_mpiRank << "\n";
  std::cout << "++++++++++++++++++++++++\n\n";


  m_fnameDmpGlbDOFs =
    "dmpGlbDOFs_" + std::to_string(m_mpiRank) + ".txt";

  m_fnameDmpLocLinOperator = 
    "dmpLocLinOperator_" + std::to_string(m_mpiRank) + ".txt";
  m_fnameDmpLocRHS =
    "dmpLocRHS_" + std::to_string(m_mpiRank) + ".txt";
  m_fnameDmpGlbDirichletInd =
    "dmpGlbDirichletIds_" + std::to_string(m_mpiRank) + ".txt";



}

Data::~Data()
{
  if (m_p_interfaceOperatorB) delete m_p_interfaceOperatorB;
}

void Data::SymbolicAssembling(std::vector<int>& elemntIds)
{
//  m_container.insert(m_container.end(), elemntIds.begin(),
 //     elemntIds.end());
  SymbolicAssembling(elemntIds.data(),elemntIds.size());
}
void Data::SymbolicAssembling(int vals[],int size)
{

#if defined(DMP_INPUTS)
  m_DumpInputsGlbDOFsPerElements(vals, size);
#endif

  // gather global DOFs of subdomain elemnts to create
  // mapping vector
  m_container.insert(m_container.end(), vals, vals + size);
}

void Data::Finalize()
{
  std::cout.rdbuf(m_p_backup);        // restore cout's original streambuf

  m_filestr.close();
}

int Data::FinalizeSymbolicAssembling()
{

#if defined(DMP_INPUTS)
  if (m_ofstrDumpInputsGlbDOFs.is_open())
    m_ofstrDumpInputsGlbDOFs.close();
#endif

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


  int dumpMatrices = m_root.get<int>("outputs.dumpMatrices",0);
  if (dumpMatrices > 0)
  {
    m_p_interfaceOperatorB->printInterfaceDOFs("B_");
    m_p_interfaceOperatorB->printNeighboursRanks("neighbRanks_");
  }


  return m_domain.GetNumberOfPrimalDOFs();

}



void Data::NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS)
{

  int ndofs = glbIds.size();
  int nK = valLocK.size();
  int nRHS = glbIds.size();
  if (nRHS != ndofs && (pow(ndofs,2) != nK))
    std::runtime_error("NumericalAssembling inputs have wrong dimensions");
  NumericAssembling(glbIds.data(),valLocK.data(),valLocRHS.data(),ndofs);

}


void Data::NumericAssembling(int glbIds[],
    double valLocK[], double valLocRHS[], int nDofs)
{

#if defined(DMP_INPUTS)
  m_DumpInputsKfPerElements(glbIds, valLocK, valLocRHS, nDofs);
#endif

  m_domain.GetStiffnessMatrix()->AddElementContribution(glbIds, valLocK, nDofs);
  m_domain.GetStiffnessMatrix()->AddRHSContribution(glbIds, valLocRHS,nDofs);
}

void Data::FinalizeNumericAssembling()
{

#if defined(DMP_INPUTS)
  if (m_ofstrDumpInputsLinOperator.is_open())
    m_ofstrDumpInputsLinOperator.close();
  if (m_ofstrDumpInputsRHS.is_open())
    m_ofstrDumpInputsRHS.close();
#endif


  // DIRICHLET BOUNDARY CONDTION
  if (m_DirichletGlbDofs.size() == 0)
    std::runtime_error("Dirichlet BC must be provided \
        before \"FinalizeNumericAssembling\" is called");

  m_domain.SetDirichletDOFs(m_DirichletGlbDofs);


  // ASSEMBLE & FACTORIZE LINEAR OPERATOR
  m_domain.GetStiffnessMatrix()->FinalizeNumericPart(m_domain.GetDirichletDOFs());

  int dumpMatrices = m_root.get<int>("outputs.dumpMatrices",0);
  if (dumpMatrices!=0){
    m_domain.GetStiffnessMatrix()->PrintStiffnessMatrix("K_");
    m_domain.GetStiffnessMatrix()->PrintKernel("R_");
    m_domain.GetStiffnessMatrix()->PrintRHS("f_");
    m_domain.GetStiffnessMatrix()->PrintNullPivots("nullPivots_");
  }

  // DIRICHLET PRECONDITIONER
  m_domain.HandlePreconditioning();

  // NUMEBRING FOR GtG MATRIX
  m_SetKernelNumbering();

  // FETI COARSE SPACE 
  auto mpB = m_p_interfaceOperatorB;
  mpB->FetiCoarseSpace(_GetDefectPerSubdomains());

}


void Data::SetDirichletDOFs(int vals[],int _size)
{
  m_DirichletGlbDofs.resize(_size);
  for (int dof = 0; dof < _size; dof++)
    m_DirichletGlbDofs[dof] = vals[dof];

#if DMP_INPUTS
  m_DumpInputsDirichletIds(vals,_size);
#endif
}

void Data::SetDirichletDOFs(std::vector<int>&v)
{
  //
  int _size = (int)v.size();
  SetDirichletDOFs(v.data(),_size);

//  m_DirichletGlbDofs.resize(_size); 
//  //
//  for (int dof = 0; dof < _size; dof++)
//    m_DirichletGlbDofs[dof] = v[dof];
}



void Data::m_SetKernelNumbering()
{


  int mpisize = m_domain.GetMpiSize();

  int myDefect = m_domain.GetStiffnessMatrix()->GetDefect();
  m_defectPerSubdomains.resize(mpisize,myDefect);


  if (m_verboseLevel>2) 
  {
    std::cout << " m_defectPerSubdomains before AlltoallInt..." << std::endl;
    for (auto& iw : m_defectPerSubdomains) std::cout<< iw << ' ';
    std::cout << '\n';
  }



  int one = 1;
  m_domain.hmpi.AlltoallInt(m_defectPerSubdomains.data(),mpisize,mpisize);




  if (m_verboseLevel>2) 
  {
    std::cout << "mpisize:  " << mpisize  << std::endl;
    std::cout << "myDefect: " << myDefect << std::endl;
    std::cout << " m_defectPerSubdomains ..." << std::endl;
    for (auto& iw : m_defectPerSubdomains) std::cout<< iw << ' ';
    std::cout << '\n';
  }

  m_container.resize(0);
  m_container.shrink_to_fit();
  DBGPRINT
}


void Data::Solve(double vals[], int nrows, int ncols)
{

  Eigen::Map<Eigen::MatrixXd> solution(vals, nrows, ncols);
  Solve(solution);

}



void Data::Solve(Eigen::Ref<Eigen::MatrixXd> solution)
{


  bool solverState = m_solver.pcpg(*this,solution);


  if (solverState)
    std::cout << "successful\n";
  else
    std::cout << "Iterative solver didn't finish successfully.\n";

  int dumpMatrices = m_root.get<int>("outputs.dumpMatrices",0);
  if (dumpMatrices!=0)
  {
    std::string fname = "solution_" + std::to_string(m_mpiRank) + ".txt";
    tools::printMatrix(solution,fname);

  }
}



void Data::_dumpTxtFiles(Eigen::Ref<Eigen::MatrixXd> solution)
{
//  m_domain.GetStiffnessMatrix()->PrintStiffnessMatrix("K_");
//  m_domain.GetStiffnessMatrix()->PrintKernel("R_");
//  m_domain.GetStiffnessMatrix()->PrintRHS("f_");
//  m_domain.GetStiffnessMatrix()->PrintNullPivots("nullPivots_");
//  m_p_interfaceOperatorB->printInterfaceDOFs("B_");
//  m_p_interfaceOperatorB->printNeighboursRanks("neighbRanks_");
//
//  std::string fname = "solution_" + std::to_string(m_mpiRank) + ".txt";
//  tools::printMatrix(solution,fname);
}





///////////////////////////////////////
// build-in assembler for fast testing
///////////////////////////////////////




void Data::howMuchAssembled(int iCell, int nel){

  double aaa = static_cast<double> (100 * (iCell + 1)) / nel;
  if ((fmod(aaa, 5)) < 1e-2)
    std::cout << aaa << "% " << "(" << iCell << "," << nel <<
      ") ... fmod(aaa,5)" << (fmod(aaa, 5)) << std::endl;

}

void Data::m_dbg_printStats()
{
//// test
//  auto mpIO = m_p_interfaceOperatorB;
//
//  int neqDual = m_domain.GetNumberOfDualDOFs();
//  int neqPrimal = m_domain.GetNumberOfPrimalDOFs();
//
//  int nRHS = 3;
//
//  Eigen::MatrixXd vecD(neqDual,nRHS);
//  for (int col = 0; col < nRHS; col++)
//    vecD.col(col) = 10.0 * (col + 1) * Eigen::VectorXd::Ones(neqDual);
//  Eigen::MatrixXd vecP(neqPrimal,nRHS);
//
//  std::cout << "dual\n";
//  for (int i = 0 ; i < neqDual; i++)
//  {
//    for (int j = 0 ; j < nRHS; j++)
//    {
//      std::cout << vecD(i,j) << ' ';
//    }
//    std::cout << std::endl;
//  }
//  std::cout << " *" << std::endl;
//
//
//
//  mpIO->multBt(vecD,vecP);
//
//  std::cout << "primal\n";
//  for (int i = 0 ; i < neqPrimal; i++)
//  {
//    for (int j = 0 ; j < nRHS; j++)
//    {
//      std::cout << vecP(i) << ' ';
//    }
//    std::cout << std::endl;
//  }
//  std::cout << " *" << std::endl;
//
//  mpIO->multB(vecP,vecD);
//  
//  std::cout << "dual\n";
//  for (int i = 0 ; i < neqDual; i++)
//  {
//    for (int j = 0 ; j < nRHS; j++)
//    {
//      std::cout << vecD(i,j) << ' ';
//    }
//    std::cout << std::endl;
//  }
//  std::cout << " *" << std::endl;


}




void Data::PathToSolverOptionFile(std::string path2file)
{

  std::cout << "parsing \"hddConf.json\" file\n";
  pt::read_json(path2file, m_root);

  m_verboseLevel = m_root.get<int>("outputs.verbose",0);

//  std::string precond = m_root.get<std::string>("solver.preconditioner");
//  std::cout << "preconditioner: " << precond << std::endl;
//  int nx = m_root.get<int>("meshgenerator.numberOfElements_x",0);
//  std::cout << "meshgenerator.numberOfElements_x: " << nx << std::endl;

}



void Data::m_DumpInputsKfPerElements(int glbIds[], double valLocK[],
    double valLocRHS[], int nDofs)
{
  // local linear operator
  if (!m_ofstrDumpInputsLinOperator.is_open())
    m_ofstrDumpInputsLinOperator.open(m_fnameDmpLocLinOperator);

  for (int dof = 0; dof < nDofs; dof++)
  {
    m_ofstrDumpInputsLinOperator << std::setprecision(16) << " " << 
      glbIds[dof];
  }
  m_ofstrDumpInputsLinOperator << std::endl;

  for (int dof = 0; dof < nDofs*nDofs; dof++)
    m_ofstrDumpInputsLinOperator << " " << valLocK[dof];
  m_ofstrDumpInputsLinOperator << std::endl;

  // local RHS
  if (!m_ofstrDumpInputsRHS.is_open())
    m_ofstrDumpInputsRHS.open(m_fnameDmpLocRHS);

  for (int dof = 0; dof < nDofs; dof++)
    m_ofstrDumpInputsRHS << std::setprecision(16) << " " << glbIds[dof];
  m_ofstrDumpInputsRHS << std::endl;

  for (int dof = 0; dof < nDofs; dof++)
    m_ofstrDumpInputsRHS << " " << valLocRHS[dof];
  m_ofstrDumpInputsRHS << std::endl;

}

void Data::m_DumpInputsDirichletIds(int vals[],int _size)
{
  if (!m_ofstrDumpInputsDirichletInd.is_open())
    m_ofstrDumpInputsDirichletInd.open(m_fnameDmpGlbDirichletInd);

  for (int dof = 0; dof < _size; dof++)
    m_ofstrDumpInputsDirichletInd << " " << vals[dof];

  m_ofstrDumpInputsDirichletInd.close();
}


void Data::m_DumpInputsGlbDOFsPerElements(int glbIds[], int nDofs)
{
  // local linear operator
  if (!m_ofstrDumpInputsGlbDOFs.is_open())
    m_ofstrDumpInputsGlbDOFs.open(m_fnameDmpGlbDOFs);

  for (int dof = 0; dof < nDofs; dof++)
    m_ofstrDumpInputsGlbDOFs << " " << glbIds[dof];
  m_ofstrDumpInputsGlbDOFs << std::endl;
}
