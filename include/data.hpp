#pragma once

#include "mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "types.hpp"




class Data{
  public:  
    Data()  {}
    ~Data(){}
    void assembly_elasticity(Mesh*);
    void setDirichlet(SpMat &mat, 
      std::vector<int> dirInd);
      SpMat& getK(){return m_K;}

  private:
    SpMat m_K;
    Eigen::VectorXd m_rhs;

};
