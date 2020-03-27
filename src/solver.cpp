#include "../include/solver.hpp"
#include <iostream>
#include "../include/types.hpp"
//#include <Eigen/PardisoSupport>


using namespace std;

Solver::Solver(Mesh* mesh,Data* data)
{



  std::vector<int> dirDOFs = mesh->getDirDOFs();

//  for (auto& id : dirDOFs) cout << id << endl;

  SpMat& K = data->getK();
//  data->svd0(K);
  data->printMatrix(K,"K0.txt"); 
  data->setDirichlet(K, dirDOFs);
//  data->svd0(K);
  data->printMatrix(K,"K1.txt"); 
  

  Eigen::VectorXd& rhs = data->getRHS();


  double sum(0);
  for (int dof = 0; dof < rhs.size(); dof++)
    sum += rhs(dof);
  cout << "sum(rhs) = " << sum - 9.81 << endl;

  
  data->setDirichlet(rhs, dirDOFs);

  m_directSolver.sym_factor(K);
  m_directSolver.num_factor();


  Eigen::VectorXd solution(rhs.size());

  m_directSolver.solve(rhs,solution);
  mesh->addSolution(solution);





  




}
