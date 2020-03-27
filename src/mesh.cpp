#include "../include/mesh.hpp"

#include <iostream>
#include <mpi.h>
#include <vtkSmartPointer.h>
#include <vector>
#include <vtkIdList.h>
#include <vtkQuadraticQuad.h> 
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>

using namespace std;


Mesh::~Mesh(){}


void Mesh::print()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
  if (m_rank == 0)
  {
    std::cout << "m_ne = " << m_ne << std::endl;
    std::cout << "m_ns = " << m_ns << std::endl;
    std::cout << "m_nl = " << m_nl << std::endl;
  }
}





int Mesh::generateMesh(int ne, int ns, int nl)
{

  int m_ne = ne;
  int m_ns = ns;
  int m_nl = nl;

  std::string fname = "test.vtu";
  m_mesh = vtkUnstructuredGrid::New();
   // vtkSmartPointer<vtkUnstructuredGrid>::New();
  double Lx(1), Ly(1);
  

  int nel = m_ne * m_ns * m_nl;
 
  m_mesh = squareMesh(Lx,Ly,nel,nel);


  return 0;

}

vtkUnstructuredGrid* Mesh::squareMesh(double Lx2D, double Ly2D, int nex2D, int ney2D)
{

  int nPoints = (nex2D + 1) * (ney2D + 1);// +ne1D + 1;
  int nCells = nex2D * ney2D;// +ne1D;

  std::vector<int> contactNodes;
  std::vector<int> DirichletNodes;



  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();


	double _x, _y;
	double dx2D = Lx2D / nex2D;
	double dy2D = Ly2D / ney2D;
	int nPoints_ctrl = 0;
	for (int j = 0; j < ney2D + 1; j++){
		for (int i = 0; i < nex2D + 1; i++){
			_x = i * dx2D;
			_y = j * dy2D;
			points->InsertNextPoint(_x, _y, 0);

			if (j == 0)
				contactNodes.push_back(nPoints_ctrl);
			if (j == ney2D)
				DirichletNodes.push_back(nPoints_ctrl);

			nPoints_ctrl++;

		}
	}


	int nnodsBasicGrid = nPoints_ctrl;

	int ni;
	for (int j = 0; j < 2 * ney2D + 1; j++){
		ni = nex2D + j % 2;
		for (int i = 0; i < ni; i++){
			_x = i * dx2D + static_cast<double>((j % 2) == 0) * dx2D * 0.5;
			_y = j * dy2D * 0.5;
			points->InsertNextPoint(_x, _y, 0);

			if (j == 0)
				contactNodes.push_back(nPoints_ctrl);
			if (j == 2 * ney2D)
				DirichletNodes.push_back(nPoints_ctrl);

			nPoints_ctrl++;
		}
	}



  int nDirNds = DirichletNodes.size();
  m_DirichletDofs.resize(2 * nDirNds);
  for (int inod = 0 ; inod < nDirNds; inod++)
  {
    m_DirichletDofs[2 * inod + 0] = 2 * DirichletNodes[inod] + 0;
    m_DirichletDofs[2 * inod + 1] = 2 * DirichletNodes[inod] + 1;
  }

//	double *xy;
//	for (int i = 0; i < points->GetNumberOfPoints(); i++)
//  {
//		xy = points->GetPoint(i);
//		points->SetPoint(i, xy);
//	}

//	vtkSmartPointer<vtkUnstructuredGrid> ug = vtkSmartPointer<vtkUnstructuredGrid>::New();

  vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
	ug->SetPoints(points);



	vtkSmartPointer<vtkIdList> quad_ids = vtkSmartPointer<vtkIdList>::New();

	quad_ids->Allocate(8, 0);

	ug->Allocate(nCells, 1);

	vector <int> tmp_vec(4, 0);
	for (int j = 0; j < ney2D; j++){
		for (int i = 0; i < nex2D; i++){
			tmp_vec[0] = i + 0;
			tmp_vec[1] = i + 1;
			tmp_vec[2] = i + 1 + (nex2D + 1);
			tmp_vec[3] = i + 0 + (nex2D + 1);
			for (int ii = 0; ii < tmp_vec.size(); ii++)
				quad_ids->InsertId(ii, tmp_vec[ii] + (nex2D + 1) * j);

			tmp_vec[0] = i + 0;
			tmp_vec[1] = i + nex2D + 1;
			tmp_vec[2] = i + 2 * nex2D + 1;
			tmp_vec[3] = i + nex2D;
			for (int ii = 0; ii < tmp_vec.size(); ii++)
				quad_ids->InsertId(ii + tmp_vec.size(), 
            tmp_vec[ii] + (2 * nex2D + 1) * j + nnodsBasicGrid);

			ug->InsertNextCell(VTK_QUADRATIC_QUAD, quad_ids);
		}
	}

	{
		vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
		intArray->SetName("FormulationId");
		intArray->SetNumberOfTuples(nCells);
		intArray->SetNumberOfComponents(1);
		intArray->Fill(1);
		ug->GetCellData()->AddArray(intArray);
	}
	{
		vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
		intArray->SetName("MaterialId");
		intArray->SetNumberOfTuples(nCells);
		intArray->SetNumberOfComponents(1);
		intArray->Fill(1);
		ug->GetCellData()->AddArray(intArray);
	}

  return ug;
}



void Mesh::writeMesh(vtkUnstructuredGrid *vtkVolumeMesh,
    std::string filename, bool asciiOrBinaryVtu){
	
	vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
		vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(filename.c_str());

	if (asciiOrBinaryVtu){
		writer->SetDataModeToAscii();
	}
	else{
		writer->SetDataModeToBinary();
	}

	writer->SetInputData(vtkVolumeMesh);
	writer->Write();
}

void Mesh::addSolution(Eigen::VectorXd& solutionXd)
{
 		vtkSmartPointer<vtkDoubleArray> 
      solutionE = vtkSmartPointer<vtkDoubleArray>::New();
		solutionE->SetName("displacement");
		solutionE->SetNumberOfComponents(3);
		solutionE->SetNumberOfTuples(m_mesh->GetNumberOfPoints());
		for (auto i = 0; i < m_mesh->GetNumberOfPoints(); i++) {
			solutionE->SetTuple3(i, 
          solutionXd(2 * i + 0), solutionXd(2 * i + 1), 0);

		}
		m_mesh->GetPointData()->AddArray(solutionE);
	}




