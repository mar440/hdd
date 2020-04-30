#include <iostream>
#include <mpi.h>
#include <vtkUnstructuredGrid.h>
#include "include/mesh.hpp"
#include "include/data.hpp"
#include "include/solver.hpp"
#include "include/element.hpp"
#include "include/linearAlgebra.hpp"
#include "include/domain.hpp"
#include "include/types.hpp"

#include <vtkCell.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>



using namespace Eigen;

int main(int argc, char **argv) 
{

  MPI_Init(&argc, &argv);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);

  int rank;

  MPI_Comm_rank(comm, &rank);
//## HDD ##  
//# initialize library (copy mpi_comm)
  Data dataH(&comm);
  dataH.ParseJsonFile("../hddConf.json");

//////////////////////
// Built-in mesher and
// assembler of A*x=b
  Mesh mesh;
  {
    // build-in generator
    auto meshOpts = dataH.GetChild("builtInTestCaseSetting.mesh");
    mesh.GenerateMesh(rank,meshOpts);
  }

// decomposition
  mesh.extractSubdomainMesh();


  // TODO make it inside the HDD library
  dataH.GetDomain()->SetNeighboursRanks(mesh.getNeighboursRanks());


  // symbolic part
  {
    vtkUnstructuredGrid* subMesh = mesh.getSubdomainMesh();
    auto glbNum = subMesh->GetCellData()->GetArray(GLOBAL_NUMBERING);
    int nElemeSubdomain = subMesh->GetNumberOfCells();
    int nP = subMesh->GetCell(0)->GetNumberOfPoints();

    std::vector<int> glbIds;
    glbIds.resize(2 * nP);

    for (int iE = 0; iE < nElemeSubdomain; iE++)
    {
      for (int iP = 0; iP < nP; iP++)
      {
        int iGI = (int) glbNum->GetComponent(iE,iP);
        glbIds[2 * iP    ] = 2 * iGI;
        glbIds[2 * iP + 1] = 2 * iGI + 1;
      }
//## HDD ##
//#  loop over subdomain elements
//#  passing global DOFs element by element - metis decomposition
      dataH.SymbolicAssembling(glbIds);
    }
  }

//## HDD ##
//# global DOFs where Dirichlet BC is applied
  dataH.SetDirichletDOFs(mesh.getDirDOFs());

//## HDD ##
//# creation mapping vectors etc ...
  dataH.FinalizeSymbolicAssembling();


  {
    vtkUnstructuredGrid* subMesh = mesh.getSubdomainMesh();
    auto glbNum = subMesh->GetCellData()->GetArray(GLOBAL_NUMBERING);
    int nElemeSubdomain = subMesh->GetNumberOfCells();
    int nP = subMesh->GetCell(0)->GetNumberOfPoints();

    std::cout << "np: " << nP << '\n';

    std::vector<int> glbIds;
    glbIds.resize(2 * nP);

    MatrixXd K_loc;
    VectorXd f_loc;


    auto matOpts = dataH.GetChild("builtInTestCaseSetting.material");

    double mat_E0     = matOpts.get<double>("YoungsModulus",0);
    double mat_mu     = matOpts.get<double>("poissonRatio",0);
    double mat_rho    = matOpts.get<double>("density",0);
    double mat_ratio  = matOpts.get<double>("ratio",0);



    auto MatId = subMesh->GetCellData()->GetArray(MATERIAL_ID);
    for (int iE = 0; iE < nElemeSubdomain; iE++)
    {
      auto cell = subMesh->GetCell(iE);
      Element *element;

      switch (cell->GetCellType()){
        case VTK_QUADRATIC_QUAD:
          element = new QUADRATIC_QUAD;
          break;
        case VTK_QUAD:
          element = new QUAD;
          break;
        default:
          continue;
      }

      int matLabel = MatId->GetTuple1(iE);

      double weight = (mat_ratio - 1.0) * matLabel + 1.0;

      double YoungModulus = mat_E0 * weight;
      double rho = mat_rho * weight;
      double mat_E_mu_rho[3] = {YoungModulus, mat_mu, rho};

      element->assembly_elasticity(K_loc, f_loc, cell, mat_E_mu_rho);

#if DBG > 4
      auto sv = linalg::svd0(K_loc);
      for (int i = 0; i < sv.size(); i++)
        std::cout << sv(i) << ' ';
      std::cout << '\n';
#endif


      for (int iP = 0; iP < nP; iP++)
      {
        int iGI = (int) glbNum->GetComponent(iE,iP);
        glbIds[iP     ] = 2 * iGI;
        glbIds[iP + nP] = 2 * iGI + 1;
      }
      std::vector<double> stdK(K_loc.data(),
          K_loc.data() + K_loc.size());
      std::vector<double> stdF(f_loc.data(),
          f_loc.data() + f_loc.size());

//## HDD ##
//# loop over elements - passing local stiffness matrices
      dataH.NumericAssembling(glbIds,stdK,stdF);
    }
  }
//## HDD ##
//# numerical factorization etc.
  dataH.FinalizeNumericAssembling();



//## HDD ##
//# solving 
  Eigen::VectorXd solution;
  dataH.Solve(solution);


  auto otpOpts = dataH.GetChild("outputs.builtInTestCase");
  bool saveBuiltInMesh = otpOpts.get<bool>("saveEachSubdomainMesh",false);
  std::cout << "saveBuiltInMes: " << saveBuiltInMesh << '\n';
  if (saveBuiltInMesh)
  {
    mesh.addSolution(mesh.getSubdomainMesh(),solution);
    mesh.SaveDecomposedMesh();
  }


//## HDD ##
//# finalizing 
  dataH.Finalize();

  MPI_Comm_free(&comm);
  MPI_Finalize();

  return 0;
}

