#pragma once
#include "directSolver.hpp"
#include <Eigen/Dense>


class Data;

class Solver
{
  public:
    Solver(){}
    ~Solver(){}
    bool pcpg(Data& data, Eigen::Ref<Eigen::MatrixXd> _solution);
    bool mpcpg(Data&, Eigen::VectorXd&);
  private:

};
