#pragma once
#include "data.hpp"
#include "mesh.hpp"
#include "directSolver.hpp"


class Solver
{
  public:
    Solver(Data* data);
    ~Solver(){}
    void pcpg();
  private:
    DirectSolver m_directSolver;

};
