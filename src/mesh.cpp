#include "../include/mesh.hpp"
#include "../include/linearAlgebra.hpp"

#include <iostream>
//#include <mpi.h>
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
#include <math.h>

#include <metis.h>
//#include "/home/mar440/usr/src/metis-5.1.0/include/metis.h"




using namespace std;


Mesh::~Mesh()
{
  if (m_mesh) m_mesh->Delete();
  if (m_subdomainMesh) m_subdomainMesh->Delete();

}

int Mesh::GenerateMesh(int _rank, boost::property_tree::ptree meshOptions)
{


  m_rank = _rank;



  m_numbOfStrips    = meshOptions.get<double>("numberOfStrips",3);

  m_nex             = meshOptions.get<int>("numberOfElements_x",0);
  m_ney             = meshOptions.get<int>("numberOfElements_y",0);
  m_nsxOneSub       = meshOptions.get<int>("numberOfSubdomains_x",0);
  m_nsyOneSub       = meshOptions.get<int>("numberOfSubdomains_y",0);
  m_Lx2D            = meshOptions.get<double>("lenght_x",0);
  m_Ly2D            = meshOptions.get<double>("lenght_y",0);
  m_x0              = meshOptions.get<double>("shift_x",0);
  m_y0              = meshOptions.get<double>("shift_y",0);
  int nlevel        = meshOptions.get<int>("numberOfLevels",0);
  std::string st0   = meshOptions.get<string>("typeOfElement");
  
  if (st0 == "linear")    m_typeOfElement = LINEAR;
  else if (st0 == "quad") m_typeOfElement = QUADRATIC;
  else
    std::runtime_error("typeOfElement may be \"linear\" or \"quad\"");






  if (m_rank == 0)
  {
    std::cout << "m_nex         " << m_nex        << std::endl;
    std::cout << "m_ney         " << m_ney        << std::endl;
    std::cout << "m_nsxOneSub   " << m_nsxOneSub  << std::endl;
    std::cout << "m_nsyOneSub   " << m_nsyOneSub  << std::endl;
    std::cout << "m_Lx2D        " << m_Lx2D       << std::endl;
    std::cout << "m_Ly2D        " << m_Ly2D       << std::endl;
    std::cout << "m_x0          " << m_x0         << std::endl;
    std::cout << "m_y0          " << m_y0         << std::endl;
  }



// multiplied by (2^levels)
  m_nsx = m_nsxOneSub * pow(2,nlevel);
  m_nsy = m_nsyOneSub * pow(2,nlevel);

  std::string fname = "test.vtu";
  m_mesh = vtkUnstructuredGrid::New();

  m_mesh = squareMesh();


  return 0;

}

vtkUnstructuredGrid* Mesh::squareMesh()
{

  int nex2D = m_nex * m_nsx;
  int ney2D = m_ney * m_nsy;

  int nCells = nex2D * ney2D;

  std::vector<int> nodesOnBottom;
  std::vector<int> nodesOnRight;

  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();





  double _x, _y;
  double dx2D = m_Lx2D / nex2D;
  double dy2D = m_Ly2D / ney2D;
  int nPoints_ctrl = 0;
  for (int j = 0; j < ney2D + 1; j++){
    for (int i = 0; i < nex2D + 1; i++){
      _x = i * dx2D + m_x0;
      _y = j * dy2D + m_y0;
      points->InsertNextPoint(_x, _y, 0);

      if (j == 0)
        nodesOnBottom.push_back(nPoints_ctrl);
      if (i == nex2D)
        nodesOnRight.push_back(nPoints_ctrl);

      nPoints_ctrl++;

    }
  }

  int nnodsBasicGrid = nPoints_ctrl;

  if (m_typeOfElement == QUADRATIC)
  {
    int ni;
    for (int j = 0; j < 2 * ney2D + 1; j++){
      ni = nex2D + j % 2;
      for (int i = 0; i < ni; i++){
        _x = i * dx2D + static_cast<double>((j % 2) == 0) * dx2D * 0.5 + m_x0;
        _y = j * dy2D * 0.5 + m_y0;
        points->InsertNextPoint(_x, _y, 0);

        if (j == 0)
          nodesOnBottom.push_back(nPoints_ctrl);
        if (i == 2 * nex2D)
          nodesOnRight.push_back(nPoints_ctrl);

        nPoints_ctrl++;
      }
    }
  }



  int nDirNdsRight  = nodesOnRight.size();
  int nDirNdsBottom = nodesOnBottom.size();
  m_DirichletDofs.resize(2 * nDirNdsRight + 2 * nDirNdsBottom);

  nPoints_ctrl = 0;

  for (int inod = 0 ; inod < nDirNdsRight; inod++)
  {
    m_DirichletDofs[2 * nPoints_ctrl + 0] = 2 * nodesOnRight[inod] + 0;
    m_DirichletDofs[2 * nPoints_ctrl + 1] = 2 * nodesOnRight[inod] + 1;
    nPoints_ctrl++;
  }
  for (int inod = 0 ; inod < nDirNdsBottom; inod++)
  {
    m_DirichletDofs[2 * nPoints_ctrl + 0] = 2 * nodesOnBottom[inod] + 0;
    m_DirichletDofs[2 * nPoints_ctrl + 1] = 2 * nodesOnBottom[inod] + 1;
    nPoints_ctrl++;
  }


  vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
  ug->SetPoints(points);
  vtkSmartPointer<vtkIdList> quad_ids = vtkSmartPointer<vtkIdList>::New();
  ug->Allocate(nCells, 1);


  if (m_typeOfElement == LINEAR)
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

              if (m_typeOfElement == LINEAR)
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

    int nel = ug->GetNumberOfCells();
    for (int ie = 0; ie < nel; ie++)
    {
      auto cell = ug->GetCell(ie);
      int np = cell->GetNumberOfPoints();
      auto cellPoints =  cell->GetPoints();
      double xyz[] = {0,0,0};
      for (int ip = 0 ; ip < np; ip++){
        auto point = cellPoints->GetPoint(ip);
        xyz[0] += point[0] / np;
        xyz[1] += point[1] / np;
        xyz[2] += point[2] / np;
      }

      double meanL = 0.5 * (m_Lx2D + m_Ly2D);
      double meanXY = xyz[1] + xyz[0];

      double pi = asin(1.) * 2; 

      double ratio = sin(pi * meanXY  * m_numbOfStrips)/ (meanL);


      if (ratio > 0) intArray->SetTuple1(ie,1);
      else intArray->SetTuple1(ie,0);
    }

  }

//  writeMesh(ug,"squareMesh.vtu",true);
  return ug;
}


void Mesh::SaveDecomposedMesh()
{

    std::string fnameVtk = "mesh_" + std::to_string(m_rank) + ".vtu";
    writeMesh(m_subdomainMesh, fnameVtk, true);
}

void Mesh::writeMesh(vtkUnstructuredGrid *vtkVolumeMesh,
    std::string filename, bool asciiOrBinaryVtu)
{

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

void Mesh::addSolution(
    vtkUnstructuredGrid* ug, Eigen::VectorXd& solutionXd)
{

  int nP = (int) ug->GetNumberOfPoints();

  if (solutionXd.size() != nP * 2)
    std::runtime_error("dimension problem");

  vtkSmartPointer<vtkDoubleArray> 
    solutionE = vtkSmartPointer<vtkDoubleArray>::New();
  solutionE->SetName(DISPLACEMENT);
  solutionE->SetNumberOfComponents(3);
  solutionE->SetNumberOfTuples(nP);
  for (auto i = 0; i < ug->GetNumberOfPoints(); i++) {
    solutionE->SetTuple3(i, 
        solutionXd(2 * i + 0), solutionXd(2 * i + 1), 0);

  }
  ug->GetPointData()->AddArray(solutionE);
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


void Mesh::ExtractSubdomainMesh()
{


  int nPmax = createMappingVectors();
  int numberOfGlobalCells = m_mesh->GetNumberOfCells();
  // ---------------------------------------------

  m_subdomainMesh = vtkUnstructuredGrid::New();
  vtkSmartPointer<vtkIdList> outputCellId = vtkSmartPointer<vtkIdList>::New();
  outputCellId->Allocate(nPmax);

  auto subId = m_mesh->GetCellData()->GetArray(SUBDOMAIN_ID);


  // count cells on subdomain
  int nSubCells(0);
  for (int iCell = 0; iCell < numberOfGlobalCells; iCell++)
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
//    cout << "array \"" << GLOBAL_NUMBERING  <<  "\" added in global mesh\n";
  }



  auto glbNumering =
    m_subdomainMesh->GetCellData()->GetArray(GLOBAL_NUMBERING);

  nSubCells = 0;
  for (int iCell = 0; iCell < numberOfGlobalCells; iCell++)
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


  {
    vtkSmartPointer<vtkIntArray> intArray = 
      vtkSmartPointer<vtkIntArray>::New();
    intArray->SetName(MATERIAL_ID);
    intArray->SetNumberOfComponents(1);
    intArray->SetNumberOfTuples(nSubCells);
    intArray->FillValue(1);

    int cnt(0);
    auto MatIdGlbMesh = m_mesh->GetCellData()->GetArray(MATERIAL_ID);
    for (int iCell = 0; iCell < numberOfGlobalCells; iCell++)
    {
      if ((int)subId->GetTuple1(iCell) == m_rank)
        intArray->SetTuple1(cnt++,MatIdGlbMesh->GetTuple1(iCell));
    }


    m_subdomainMesh->GetCellData()->AddArray(intArray);
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

//  std::string fname  = "dom_";
//  fname += std::to_string(m_rank);
//  fname += ".vtu";
//  writeMesh(m_subdomainMesh,fname,true);
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


void Mesh::DomainDecomposition()
{


    int _nparts = m_nsx * m_nsy;

    int nCells= m_mesh->GetNumberOfCells();
    int nPointsLocal = m_mesh->GetNumberOfPoints();

    int nPi(0);
    std::vector<int> eptr;
    eptr.push_back(0);

//
    int cnt = 0;
    for (int iCell = 0; iCell < nCells; iCell++)
    {
      vtkCell *oneCell=
        dynamic_cast<vtkCell*>(m_mesh->GetCell(iCell));
      nPi = oneCell->GetNumberOfPoints();
      for (int ii = 0; ii < nPi; ii++)
      {
        cnt += nPi;
        eptr.push_back(cnt);
      }
    }


    int* eind = new int[eptr[nCells]];
    cnt = 0;
    int CellDim = 2;

    for (int iCell = 0; iCell < nCells; iCell++)
    {
      vtkCell *oneCell=
        dynamic_cast<vtkCell*>(m_mesh->GetCell(iCell));
      nPi = oneCell->GetNumberOfPoints();

      auto innerCellPointIds = oneCell->GetPointIds();
      for (int ii = 0; ii < nPi; ii++)
      {
          eind[cnt] = innerCellPointIds->GetId(ii);
          cnt ++;
      }

    }

    int ncommon = 2; /* ncommon = 2 for all other types ... */
    if (CellDim == 2){
        ncommon = 2;
    }
    else if (CellDim == 3){
        ncommon = 3;
    }


    int options[METIS_NOPTIONS];
    cout << " ----------   METIS_NOPTIONS " << METIS_NOPTIONS << endl;
    options[METIS_OPTION_PTYPE    ] = METIS_PTYPE_RB;    // multilevel recursive bisectioning
    options[METIS_OPTION_OBJTYPE  ] = METIS_OBJTYPE_CUT; // edge-cut minimization
    options[METIS_OPTION_CTYPE    ] = METIS_CTYPE_RM;    // random matching
    options[METIS_OPTION_IPTYPE   ] = METIS_IPTYPE_GROW; // grows a bisction using a greedy strategy
    options[METIS_OPTION_RTYPE    ] = METIS_RTYPE_FM;    // FM-based cut refinement
    options[METIS_OPTION_NCUTS    ] = 1 ;//
    options[METIS_OPTION_NITER    ] = 10;/* Default value */
    options[METIS_OPTION_SEED     ] = -1;/* Seed of random algo */
    options[METIS_OPTION_UFACTOR  ] = 1;
    options[METIS_OPTION_NUMBERING] = 0; // C-style numbering
    options[METIS_OPTION_DBGLVL   ] = METIS_DBG_INFO;
    options[METIS_OPTION_CONTIG   ] = 1;

    int nparts = 1;
    int objval;


    if (_nparts > 0 ){
       nparts = _nparts;
    }

    int *epart = new int [nCells];
    int *npart = new int [eptr[nCells]];


    if (nparts > 1){
        METIS_PartMeshDual(&nCells,    // number of elements in the mesh       Y
                           &nPointsLocal,   //                                      Y
                           eptr.data(),     //                                      Y
                           eind,            //                                      Y
                           (int *)NULL,     // vwgt                                 Y
                           (int *)NULL,     // vsize                                Y
                           &ncommon,        //                                      Y
                           &nparts,         //                                      N
                           (real_t*)NULL,   // tpwgts                               Y
                           options,         //                                      Y
                           &objval,         //                                      Y
                           epart,           //                                      N
                           npart);          //                                      Y
    }
    else {
        for (int i = 0; i < nCells; i++)
            epart[i] = 0;

    }




  auto subId = m_mesh->GetCellData()->GetArray(SUBDOMAIN_ID);


  // count cells on subdomain
  for (int iCell = 0; iCell < nCells; iCell++)
    subId->SetTuple1(iCell,epart[iCell]);

  writeMesh(m_mesh,"squareMesh_metis.vtu",true);

  delete [] epart;
  delete [] npart;
  delete [] eind;

}











