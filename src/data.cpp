#include "../include/data.hpp"
#include "../include/element.hpp"
#include <iostream>
#include <vtkCell.h>
#include <string>



using namespace std;
using namespace Eigen;
typedef Triplet<double> T;



void Data::assembly_elasticity(Mesh *mesh)
{

	int nel = mesh->m_mesh->GetNumberOfCells();
  cout << "nel = " << nel << endl;
	int dim = 2;
	//double mat_E_mu[3] = {2.1e5, 0.3, 0.78500e9};
	double mat_E_mu[3] = {1.0, 0.3, 1.0};

	cout << "K assembly elasticity ..." << endl;
	int nDOFs(0);
  nDOFs = mesh->m_mesh->GetNumberOfPoints() * dim;
  m_rhs.resize(nDOFs);
	m_rhs.setZero();
  std::vector<T> trK;
  MatrixXd K_loc;
	VectorXd f_loc;
  
	for (int iCell = 0; iCell < nel; iCell++)
  {
    howMuchAssembled(iCell,nel);

		auto cell = mesh->m_mesh->GetCell(iCell);
		int nP = cell->GetNumberOfPoints();

		Element *element;

		switch (cell->GetCellType()){
		case VTK_QUADRATIC_QUAD:
			element = new QUADRATIC_QUAD;
			break;
		default:
			continue;
		}


		element->assembly_elasticity(K_loc, f_loc, cell, mat_E_mu);

//    auto singVals = svd0(K_loc);
//    for (int vls = 0; vls < singVals.size(); vls++)
//      cout << singVals(vls) << ' ';
//    cout << endl;
		//TODO the numbering for tetra 10 is (probably) reversed ...

		vtkIdList *ids = cell->GetPointIds();

 
		int i_glb, j_glb, i_loc, j_loc;
		for (int iDim = 0; iDim < dim; iDim++){
			for (int jDim = 0; jDim < dim; jDim++){
				for (int i = 0; i < nP; i++){
					i_glb = dim * ids->GetId(i) + iDim;
					i_loc = i_glb;//mesh->g2l_volume[i_glb];
					if (jDim == 0)
						m_rhs(i_loc) += f_loc(i + iDim * nP);
					for (int j = 0; j < nP; j++){
						j_glb = dim * ids->GetId(j) + jDim;
						j_loc = j_glb;// mesh->g2l_volume[j_glb];
//						trp_J22_renum.push_back(T(i_loc, j_loc, K_loc(i + iDim * nP, j + jDim * nP)));
              trK.push_back(T(i_glb, j_glb, K_loc(i + iDim * nP, j + jDim * nP)));
					}
				}
			}
		}
		delete element;
	}


	// global stiffness matrix for elasticity - initialization
	m_K.resize(nDOFs, nDOFs);
	m_K.setFromTriplets(trK.begin(), trK.end());



	cout << " elasticity matrix assembled ... " << endl;
}

void Data::setDirichlet(Eigen::VectorXd& v,
      std::vector<int> dirInd)
{

  for (int iD = 0; iD < static_cast<int>(dirInd.size()); iD++)
    v[dirInd[iD]] = 0;

}


void Data::setDirichlet(
    SpMat &mat,
    std::vector<int> dirInd)
{

	VectorXd diagMat = mat.diagonal();
	double meanVal = 1;// diagMat.sum() / diagMat.size();

	for (int iD = 0; iD < static_cast<int>(dirInd.size()); iD++){

		int kD = dirInd[iD];
		for (SpMat::InnerIterator it(mat, kD); it; ++it)
		{
			if (it.row() == it.col())
				it.valueRef() = meanVal;
			else
				it.valueRef() = 0;
		}
	}
}




void Data::printMatrix(SpMat& mat,std::string fname){
	{ 
		ofstream myfile(fname);
		if (myfile.is_open())
		{
			myfile << setprecision(16);
			for (int indJ = 0; indJ < mat.rows(); indJ++) {
				for (SpMat::InnerIterator it(mat, indJ); it; ++it) {
					myfile << it.row() << " " << it.col() << " " << it.value() << "\n";
				}
			}
			myfile.close();
		}
	}

}


VectorXd Data::svd0(Eigen::MatrixXd& dmat)
{
  JacobiSVD<MatrixXd> svd(dmat, ComputeThinU | ComputeThinV);
//  cout << "Its singular values are:" << endl << svd.singularValues() << endl;
//  cout << "Its left singular vectors are the columns of the thin U matrix:" << 
//    endl << svd.matrixU() << endl;
//  cout << "Its right singular vectors are the columns of the thin V matrix:" << 
//    endl << svd.matrixV() << endl;
//  Vector3f rhs(1, 0, 0);
//  cout << "Now consider this rhs vector:" << endl << rhs << endl;
//  cout << "A least-squares solution of m*x = rhs is:" << endl << svd.solve(rhs) << endl;
  Eigen::MatrixXd singVals = svd.singularValues();
  return singVals;
}

VectorXd Data::svd0(SpMat& smat)
{
  Eigen::MatrixXd dmat = MatrixXd(smat);
  return svd0(dmat);
}



void Data::howMuchAssembled(int iCell, int nel){

		double aaa = static_cast<double> (100 * (iCell + 1)) / nel;
		if ((fmod(aaa, 5)) < 1e-2)
			cout << aaa << "% " << "(" << iCell << "," << nel <<
        ") ... fmod(aaa,5)" << (fmod(aaa, 5)) << endl;

}
