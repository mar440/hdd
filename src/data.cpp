#include "../include/data.hpp"
#include "../include/element.hpp"
#include "../include/linearAlgebra.hpp"
#include <iostream>
#include <vtkCell.h>
#include <string>



using namespace std;
using namespace Eigen;



Data::Data(MPI_Comm* _pcomm): m_domain(_pcomm)
{
  m_pcomm = _pcomm;

  m_container.resize(0);
}



void Data::SymbolicAssembling(std::vector<int>& elemntIds)
{
  // gather global DOFs of subdomain elemnts to create
  // mapping vector
  m_container.insert(m_container.end(), elemntIds.begin(),
      elemntIds.end());
}


void Data::FinalizeSymbolicAssembling()
{
  // maping DOF: local to global
  linalg::unique(m_container);
  m_domain.SetMappingLoc2Glb(m_container);
  m_container.clear();
  m_container.shrink_to_fit();

  //TODO assemble matrix B
  m_domain.SetInterfaces();

}

void Data::NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS)
{
  m_domain.NumericAssemblingStiffnessAndRhs(glbIds, valLocK, valLocRHS);
}


void Data::FinalizeNumericAssembling()
{

  m_domain.FinalizeStiffnessMatrixAndRhs();

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
