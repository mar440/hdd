#pragma once
#include "mesh.hpp"
#include "directSolver.hpp"


class Data;

class Solver
{
  public:
    Solver();
    ~Solver(){}
    void pcpg(Data&);
  private:
//    DirectSolver m_directSolver;
//    Data* m_p_data;

};
