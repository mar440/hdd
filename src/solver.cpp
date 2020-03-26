#include "../include/solver.hpp"
#include <iostream>
#include "../include/types.hpp"
#include <Eigen/PardisoSupport>


using namespace std;

Solver::Solver(Mesh* mesh,Data* data)
{



  std::vector<int> dirDOFs = mesh->getDirDOFs();

  for (auto& id : dirDOFs)
  {
    cout << id << endl;
  }

  SpMat& K = data->getK();
  data->setDirichlet(K, dirDOFs);




  




}
