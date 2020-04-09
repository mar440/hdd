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

int main(int argc, char **argv) 
{

  MPI_Init(&argc, &argv);
  MPI_Comm comm;
  MPI_Comm_dup(MPI_COMM_WORLD,&comm);

  int rank;

  MPI_Comm_rank(comm, &rank);

  std::string fname = "out_" + std::to_string(rank) + ".txt";

  std::ofstream out(fname);
  std::streambuf *coutbuf = std::cout.rdbuf(); //save old buf
  std::cout.rdbuf(out.rdbuf());
  cout << "SUBDOMAIN id. " << rank << "\n";

//////////////////////
// Build-in mesher and
// assembler of A*x=b
  Mesh mesh;
  {
    // global mesh
    int n_elements = 10;
    int n_subdomains = 8;
    int n_levels = 0;
    // build-in generator
    mesh.generateMesh( n_elements, n_subdomains, n_levels);
  }
  cout << "after mesh \n";

// decomposition
  {
    mesh.extractSubdomainMesh();
    std::string fnameVtk = "mesh_" + std::to_string(rank) + ".vtu";
    mesh.writeMesh(mesh.getSubdomainMesh(), fnameVtk, true);
  }

//## HDD ##  
//# initialize library (copy mpi_comm)
  Data dataH(&comm);
  cout << "after dataH\n";

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
          std::cout << "quadratic_quad\n";
          break;
        case VTK_QUAD:
          element = new QUAD;
          std::cout << "quad\n";
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
//# loop over elements - passing local stiffness matrices
      dataH.NumericAssembling(glbIds,stdK,stdF);
    }
  }
//## HDD ##
//# numerical factorization etc.
  dataH.FinalizeNumericAssembling();


  Solver solver;

  solver.pcpg(dataH);

//  std::string fnameVtk = "mesh_" + std::to_string(rank) + ".vtu";
//  mesh.writeMesh(mesh.getSubdomainMesh(), fnameVtk, true);
//
//  if (rank == 0)
//  {
//    fnameVtk = "globalMesh.vtu";
//    mesh.writeMesh(mesh.getGlobalMesh(), fnameVtk, true);
//  }

  std::cout.rdbuf(coutbuf); //reset to standard output again

  MPI_Comm_free(&comm);

  MPI_Finalize();

  return 0;
}

