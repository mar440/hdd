#pragma once
#include "data.hpp"
#include "mesh.hpp"
#include "directSolver.hpp"


class Solver
{
  public:
    Solver(Mesh* mesh,Data* data);
    ~Solver(){}
  private:
    DirectSolver m_directSolver;

};
