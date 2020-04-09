#include "../include/mesh.hpp"
#include "../include/linearAlgebra.hpp"

#include <iostream>
#include <mpi.h>
#include <vtkSmartPointer.h>
#include <vector>
#include <vtkIdList.h>
#include <vtkQuadraticQuad.h> 
#include <vtkQuad.h> 
#include <vtkIntArray.h>
#include <vtkCellData.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkDoubleArray.h>
#include <vtkPoints.h>
#include <vtkPointData.h>



enum TYPE_OF_ELEMENTS {LINEAR, QUADRATIC};

TYPE_OF_ELEMENTS TypeOfElements = LINEAR;
//TYPE_OF_ELEMENTS TypeOfElements = QUADRATIC;

using namespace std;

Mesh::Mesh()
{
  MPI_Comm_rank(MPI_COMM_WORLD, &m_rank);
}

Mesh::~Mesh()
{
  if (m_mesh)
    m_mesh->Delete();

  if (m_subdomainMesh)
    m_subdomainMesh->Delete();

}




int Mesh::generateMesh(int ne, int ns, int nl)
{

// number of elements in x, y direction per subdomain
  m_nex = m_ney = ne;

// number of subdomains in x, y direction per subdomain
// multiplied by (2^levels)
  m_nsx = m_nsy = ns * pow(2,nl);

  std::string fname = "test.vtu";
  m_mesh = vtkUnstructuredGrid::New();
  // vtkSmartPointer<vtkUnstructuredGrid>::New();
  double Lx(2), Ly(2), x0(-1), y0(-1);

  m_mesh = squareMesh(Lx,Ly, x0, y0);


  return 0;

}

vtkUnstructuredGrid* Mesh::squareMesh(double Lx2D, double Ly2D, double x0, double y0)
{





  int nex2D = m_nex * m_nsx;
  int ney2D = m_ney * m_nsy;

  int nCells = nex2D * ney2D;

  std::vector<int> nodesOnBottom;
  std::vector<int> nodesOnTop;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();





  double _x, _y;
  double dx2D = Lx2D / nex2D;
  double dy2D = Ly2D / ney2D;
  int nPoints_ctrl = 0;
  for (int j = 0; j < ney2D + 1; j++){
    for (int i = 0; i < nex2D + 1; i++){
      _x = i * dx2D + x0;
      _y = j * dy2D + y0;
      points->InsertNextPoint(_x, _y, 0);

      if (j == 0)
        nodesOnBottom.push_back(nPoints_ctrl);
      if (j == ney2D)
        nodesOnTop.push_back(nPoints_ctrl);

      nPoints_ctrl++;

    }
  }

  int nnodsBasicGrid = nPoints_ctrl;

  if (TypeOfElements == QUADRATIC)
  {
    int ni;
    for (int j = 0; j < 2 * ney2D + 1; j++){
      ni = nex2D + j % 2;
      for (int i = 0; i < ni; i++){
        _x = i * dx2D + static_cast<double>((j % 2) == 0) * dx2D * 0.5 + x0;
        _y = j * dy2D * 0.5 + y0;
        points->InsertNextPoint(_x, _y, 0);

        if (j == 0)
          nodesOnBottom.push_back(nPoints_ctrl);
        if (j == 2 * ney2D)
          nodesOnTop.push_back(nPoints_ctrl);

        nPoints_ctrl++;
      }
    }
  }



  int nDirNds = nodesOnTop.size();
  m_DirichletDofs.resize(2 * nDirNds);
  for (int inod = 0 ; inod < nDirNds; inod++)
  {
    m_DirichletDofs[2 * inod + 0] = 2 * nodesOnTop[inod] + 0;
    m_DirichletDofs[2 * inod + 1] = 2 * nodesOnTop[inod] + 1;
  }


  vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
  ug->SetPoints(points);
  vtkSmartPointer<vtkIdList> quad_ids = vtkSmartPointer<vtkIdList>::New();
  ug->Allocate(nCells, 1);


  if (TypeOfElements == LINEAR)
    quad_ids->Allocate(4, 0);
  else
    quad_ids->Allocate(8, 0);

  vector <int> tmp_vec(4, 0);

  int cnt(0);
  {
    vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
    intArray->SetName(SUBDOMAIN_ID);
    intArray->SetNumberOfComponents(1);
    intArray->SetNumberOfTuples(nCells);
    intArray->FillValue(1);
    ug->GetCellData()->AddArray(intArray);
  }


  auto subId =ug->GetCellData()->GetArray(SUBDOMAIN_ID);
  int curSubId(0);

//  for (int ly = 0; ly < m_nly; ly++){
//    for (int lx = 0; lx < m_nlx; lx++){
      for (int sy = 0; sy < m_nsy; sy++){
        for (int sx = 0; sx < m_nsx; sx++){
          for (int ey = sy * m_ney ; ey < (sy + 1) * m_ney; ey++){
            for (int ex = sx * m_nex ; ex < (sx + 1) * m_nex; ex++){


              tmp_vec[0] = ex + 0;
              tmp_vec[1] = ex + 1;
              tmp_vec[2] = ex + 1 + (nex2D + 1);
              tmp_vec[3] = ex + 0 + (nex2D + 1);
              for (int ii = 0; ii < (int)tmp_vec.size(); ii++)
                quad_ids->InsertId(ii, tmp_vec[ii] + (nex2D + 1) * ey);

              if (TypeOfElements == LINEAR)
              {
                ug->InsertNextCell(VTK_QUAD, quad_ids);
              }
              else
              {
                tmp_vec[0] = ex + 0;
                tmp_vec[1] = ex + nex2D + 1;
                tmp_vec[2] = ex + 2 * nex2D + 1;
                tmp_vec[3] = ex + nex2D;
                for (int ii = 0; ii < (int) tmp_vec.size(); ii++)
                  quad_ids->InsertId(ii + tmp_vec.size(),
                      tmp_vec[ii] + (2 * nex2D + 1) * ey + nnodsBasicGrid);

                ug->InsertNextCell(VTK_QUADRATIC_QUAD, quad_ids);
              }

              subId->SetTuple1(cnt,(double)curSubId);
              cnt++;
            }
          }
          curSubId++;
        }
      }
//    }
//  }

  {
    vtkSmartPointer<vtkIntArray> intArray =
      vtkSmartPointer<vtkIntArray>::New();
    intArray->SetName(FORMULATION_ID);
    intArray->SetNumberOfComponents(1);
    intArray->SetNumberOfTuples(nCells);
    intArray->FillValue(1);
    ug->GetCellData()->AddArray(intArray);
  }

  {
    vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
    intArray->SetName(MATERIAL_ID);
    intArray->SetNumberOfComponents(1);
    intArray->SetNumberOfTuples(nCells);
    intArray->FillValue(1);
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
  solutionE->SetName(DISPLACEMENT);
  solutionE->SetNumberOfComponents(3);
  solutionE->SetNumberOfTuples(m_mesh->GetNumberOfPoints());
  for (auto i = 0; i < m_mesh->GetNumberOfPoints(); i++) {
    solutionE->SetTuple3(i, 
        solutionXd(2 * i + 0), solutionXd(2 * i + 1), 0);

  }
  m_mesh->GetPointData()->AddArray(solutionE);
}




void Mesh::print()
{
  if (m_rank == 0)
  {
    std::cout << "m_nex = " << m_nex << std::endl;
    std::cout << "m_ney = " << m_ney << std::endl;

    std::cout << "m_nsx = " << m_nsx << std::endl;
    std::cout << "m_nsy = " << m_nsy << std::endl;

    std::cout << "m_nlx = " << m_nlx << std::endl;
    std::cout << "m_nly = " << m_nly << std::endl;
  }
}

int Mesh::createMappingVectors()
{

  int nel = m_mesh->GetNumberOfCells();

  auto subId = m_mesh->GetCellData()->GetArray(SUBDOMAIN_ID);

  int nPmax(0);
  for (int iCell = 0; iCell < nel; iCell++)
  {
//    cout << "subId = " <<   (int)subId->GetTuple1(iCell)  << endl; 
//    cout << "rank  = " <<   m_rank << endl; 

    if ((int)subId->GetTuple1(iCell) == m_rank)
    {

      vtkCell *cell = dynamic_cast<vtkCell*>(m_mesh->GetCell(iCell));
      int nP = cell->GetNumberOfPoints();
      auto innerCellPointIds = cell->GetPointIds();
      for (int ii = 0; ii < nP; ii++)
        m_l2g.push_back(innerCellPointIds->GetId(ii));

      if (nP > nPmax) nPmax = nP;

    }
  }

  linalg::unique(m_l2g);

  for (int i = 0; i < (int)m_l2g.size(); i++)
    m_g2l[m_l2g[i]] = i;


  int maxNumberOfPointOnCell(nPmax);

  return maxNumberOfPointOnCell;
}


void Mesh::extractSubdomainMesh()
{


  int nPmax = createMappingVectors();
  int nel = m_mesh->GetNumberOfCells();
  // ---------------------------------------------

  m_subdomainMesh = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkIdList> outputCellId = vtkSmartPointer<vtkIdList>::New();
  outputCellId->Allocate(nPmax);

  auto subId = m_mesh->GetCellData()->GetArray(SUBDOMAIN_ID);


  // count cells on subdomain
  int nSubCells(0);
  for (int iCell = 0; iCell < nel; iCell++)
    if ((int)subId->GetTuple1(iCell) == m_rank)
      nSubCells++;


  {
    vtkCell *cell0 =
      dynamic_cast<vtkCell*>(m_mesh->GetCell(0));

    int nP = cell0->GetNumberOfPoints();
    vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
    intArray->SetName(GLOBAL_NUMBERING);
    intArray->SetNumberOfComponents(nP);
    intArray->SetNumberOfTuples(nSubCells);
    intArray->FillValue(-1);
    m_subdomainMesh->GetCellData()->AddArray(intArray);
    cout << "array \"" << GLOBAL_NUMBERING  <<  "\" added in global mesh\n";
  }



  auto glbNumering =
    m_subdomainMesh->GetCellData()->GetArray(GLOBAL_NUMBERING);

  nSubCells = 0;
  for (int iCell = 0; iCell < nel; iCell++)
  {
    if ((int)subId->GetTuple1(iCell) == m_rank)
    {
      vtkCell *cellOrigNumb =
        dynamic_cast<vtkCell*>(m_mesh->GetCell(iCell));
      int cellType = cellOrigNumb->GetCellType();
      int nP = cellOrigNumb->GetNumberOfPoints();
      auto inputCellId = cellOrigNumb->GetPointIds();
      for (int ii = 0; ii < nP; ii++)
      {
        outputCellId->InsertId(ii,m_g2l[inputCellId->GetId(ii)]);
        glbNumering->SetComponent(nSubCells,ii,inputCellId->GetId(ii));
      }

      m_subdomainMesh->InsertNextCell(cellType, outputCellId);
      nSubCells++;
    }
  }

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

  int nPointsSubdomain = static_cast<int>(m_l2g.size());
  points->Allocate(nPointsSubdomain);

  double *_x;
  for (int ii = 0; ii < nPointsSubdomain; ii++){
    _x = m_mesh->GetPoint(m_l2g[ii]);
    points->InsertNextPoint(_x);
  }

  m_subdomainMesh->SetPoints(points);

}



std::vector<int> Mesh::getNeighboursRanks()
{

  int J = m_rank / m_nsx;
  int I = m_rank - J * m_nsx;

  std::vector<int> neighboursRanks(0);
#if DBG > 0
  cout << "subdom. coords: " << I << ", " << J << '\n';
#endif

  int leftBottomRank = m_rank - 1 - m_nsx;
  int bottomRank = m_rank - m_nsx;
  int rightBottomRank = m_rank + 1 - m_nsx;
  int leftRank = m_rank - 1;
  int rightRank = m_rank + 1;
  int leftTopRank = m_rank - 1 + m_nsx;
  int topRank = m_rank + m_nsx;
  int rightTopRank = m_rank + 1 + m_nsx;


  if (I > 0)
    neighboursRanks.push_back(leftRank);
  if (I < m_nsx - 1)
    neighboursRanks.push_back(rightRank);
  if (J > 0)
    neighboursRanks.push_back(bottomRank);
  if (J < m_nsy - 1)
    neighboursRanks.push_back(topRank);
  if (I > 0 && J > 0)
    neighboursRanks.push_back(leftBottomRank);
  if (I < m_nsx - 1 && J > 0)
    neighboursRanks.push_back(rightBottomRank);
  if (I > 0 && J < m_nsy - 1)
    neighboursRanks.push_back(leftTopRank);
  if (I < m_nsx - 1 && J < m_nsy - 1)
    neighboursRanks.push_back(rightTopRank);

  std::sort(neighboursRanks.begin(),neighboursRanks.end());


  return neighboursRanks;

}
