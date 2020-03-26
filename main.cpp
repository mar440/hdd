#include <iostream>
#include <mpi.h>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <map>
#include <vtkUnstructuredGrid.h>
#include <Eigen/Dense>
#include <metis.h>
#include "include/mesh.hpp"
#include "include/data.hpp"
#include "include/solver.hpp"


using namespace std;
using namespace Eigen;

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MatrixXd m0 = MatrixXd::Zero(3,3);

  Mesh mesh;

  {
    int n_elements = 3;
    int n_subdomains = 1;
    int n_levels = 1;
    mesh.generateMesh(
        n_elements,
        n_subdomains,
        n_levels);
  }


  if (rank == 0)
	  mesh.writeMesh(mesh.m_mesh, "test.vtu", true);


  Data data;
  data.assembly_elasticity(&mesh);

  Solver solver(&mesh,&data);

  MPI_Finalize();

  return 0;
}

