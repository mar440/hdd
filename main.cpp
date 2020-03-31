#include <iostream>
#include <fstream>
#include <mpi.h>
#include <vtkUnstructuredGrid.h>
#include "include/mesh.hpp"
#include "include/data.hpp"
#include "include/solver.hpp"
#include "include/element.hpp"
#include "include/linearAlgebra.hpp"
#include "include/domain.hpp"

#include <vtkCell.h>
#include <vtkIntArray.h>
#include <vtkCellData.h>


using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD, &comm);

  int rank;
  MPI_Comm_rank(comm, &rank);

  std::string fname = "out_" + std::to_string(rank) + ".txt";

  std::ofstream out(fname);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf());
  cout << "SUBDOMAIN id. " << rank << "\n";

/////////////////////
// Own mesh builder
// and assembler of A*x=b
  Mesh mesh;
  {
    // global mesh
    int n_elements = 2;
    int n_subdomains = 2;
    int n_levels = 0;
    // build-in generator
    mesh.generateMesh( n_elements, n_subdomains, n_levels);
  }
  cout << "after mesh \n";

// decomposition
  {
    mesh.extractSubdomainMesh();
  }

  Data dataH(&comm);
  cout << "after dataH\n";

  dataH.getDomain().SetNeighboursRanks(mesh.getNeighboursRanks());

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
      dataH.SymbolicAssembling(glbIds);
    }
  }

//## HDD ##
  dataH.FinalizeSymbolicAssembling();
  cout << "after finalizeGlobalIds\n";


  {
    vtkUnstructuredGrid* subMesh = mesh.getSubdomainMesh();
    auto glbNum = subMesh->GetCellData()->GetArray(GLOBAL_NUMBERING);
    int nElemeSubdomain = subMesh->GetNumberOfCells();
    int nP = subMesh->GetCell(0)->GetNumberOfPoints();
    std::vector<int> glbIds;
    glbIds.resize(2 * nP);

    MatrixXd K_loc;
    VectorXd f_loc;
    double mat_E_mu_rho[3] = {1.0, 0.3, 1.0};

    for (int iE = 0; iE < nElemeSubdomain; iE++)
    {
      auto cell = subMesh->GetCell(iE);
      Element *element;

      switch (cell->GetCellType()){
        case VTK_QUADRATIC_QUAD:
          element = new QUADRATIC_QUAD;
          break;
        default:
          continue;
      }

      element->assembly_elasticity(K_loc, f_loc, cell, mat_E_mu_rho);

#if DBG > 4
      auto sv = linalg::svd0(K_loc);
      for (int i = 0; i < sv.size(); i++)
        cout << sv(i) << ' ';
      cout << '\n';
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
      dataH.NumericAssembling(glbIds,stdK,stdF);
    }
  }
//## HDD ##
  dataH.FinalizeNumericAssembling();


  //  Solver solver(&mesh,&dataH);
  std::string fnameVtk = "mesh_" + std::to_string(rank) + ".vtu";
  mesh.writeMesh(mesh.getSubdomainMesh(), fnameVtk, true);

  if (rank == 0)
  {
    fnameVtk = "globalMesh.vtu";
    mesh.writeMesh(mesh.getGlobalMesh(), fnameVtk, true);
  }
  
  std::cout.rdbuf(coutbuf); //reset to standard output again

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}

