#pragma once
#include "directSolver.hpp"
#include <Eigen/Dense>


class Data;

class Solver
{
  public:
    Solver();
    ~Solver(){}
    void pcpg(Data&, Eigen::VectorXd&);
  private:
//    DirectSolver m_directSolver;
//    Data* m_p_data;

};
