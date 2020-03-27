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

    void setDirichlet(Eigen::VectorXd& v,
      std::vector<int> dirInd);

    void setDirichlet(SpMat &mat, 
      std::vector<int> dirInd);
   
    SpMat& getK(){return m_K;}
    Eigen::VectorXd& getRHS(){return m_rhs;}
    void printMatrix(SpMat& mat,std::string fname);
    
    Eigen::VectorXd svd0(SpMat& mat);
    Eigen::VectorXd svd0(Eigen::MatrixXd& mat);

  private:
    SpMat m_K;
    Eigen::VectorXd m_rhs;
    void howMuchAssembled(int,int);

};
