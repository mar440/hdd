#include "../include/data.hpp"
#include "../include/element.hpp"
#include <iostream>
#include <vtkCell.h>



using namespace std;
using namespace Eigen;
typedef Triplet<double> T;



void Data::assembly_elasticity(Mesh *mesh)
{

	int nel = mesh->m_mesh->GetNumberOfCells();
  cout << "nel = " << nel << endl;
	int dim = 2;
	double mat_E_mu[2] = {2.1e3, 0.3};

	//vector <T> trp_J22_renum;
	cout << "K assembly elasticity ..." << endl;
	int maxI = 0;
	//VectorXd solution, rhsEig, delta_solution;
	int nDOFs(0);
	int dime = 3;
	//if (mesh->l2g_volume.size() == 0){

	//	auto cell = mesh->vtkVolumeMesh->GetCell(0);
	//	vtkIdList *ids = cell->GetPointIds();
	//	mesh->l2g_volume.reserve(nel * ids->GetNumberOfIds() * dim);

	//	mesh->l2gNds_volume.reserve(nel * ids->GetNumberOfIds());



	//	for (int iCell = 0; iCell < nel; iCell++){
	//		auto cell = mesh->vtkVolumeMesh->GetCell(iCell);
	//		vtkIdList *ids = cell->GetPointIds();

	//		for (int i = 0; i < ids->GetNumberOfIds(); i++){
	//			mesh->l2gNds_volume.push_back(ids->GetId(i));
	//			for (int d = 0; d < dim; d++)
	//				mesh->l2g_volume.push_back(dim * ids->GetId(i) + d);
	//		}
	//	}

	//	std::sort(mesh->l2g_volume.begin(), mesh->l2g_volume.end());
	//	std::vector<int>::iterator it;
	//	it = std::unique(mesh->l2g_volume.begin(), mesh->l2g_volume.end());
	//	mesh->l2g_volume.resize(std::distance(mesh->l2g_volume.begin(), it));

	//	std::sort(mesh->l2gNds_volume.begin(), mesh->l2gNds_volume.end());
	//	it = std::unique(mesh->l2gNds_volume.begin(), mesh->l2gNds_volume.end());
	//	mesh->l2gNds_volume.resize(std::distance(mesh->l2gNds_volume.begin(), it));




	//	for (int i = 0; i < mesh->l2g_volume.size(); i++)
	//		mesh->g2l_volume[mesh->l2g_volume[i]] = i;
	//}

//	nDOFs = static_cast<int>(mesh->l2g_volume.size());
  nDOFs = mesh->m_mesh->GetNumberOfPoints() * dim;
  m_rhs.resize(nDOFs);
	m_rhs.setZero();
  std::vector<T> trK;
  MatrixXd K_loc;
	VectorXd f_loc;
  
	for (int iCell = 0; iCell < nel; iCell++)
  {

		double aaa = static_cast<double> (100 * (iCell + 1)) / nel;
		if ((fmod(aaa, 5)) < 1e-2)
			cout << aaa << "% " << "(" << iCell << "," << nel << ") ... fmod(aaa,5)" << (fmod(aaa, 5)) << endl;

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


