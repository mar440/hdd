#pragma once

#include "mesh.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include "types.hpp"

#include "domain.hpp"
#include <mpi.h>



class Data{
  public:
    Data(MPI_Comm*);
    ~Data(){}

    void SymbolicAssembling(std::vector<int>& glbIds);
    void FinalizeSymbolicAssembling();

    void NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS);
    void FinalizeNumericAssembling();

    Domain& getDomain(){return  m_domain;}



  private:
    MPI_Comm* m_pcomm;
    std::vector<int> m_container;
    Domain m_domain;

    int m_cnt_setLocalMatrix;





///////////////////////////////////////
// build-in assembler for fast testing
///////////////////////////////////////
  public:
//    void assembly_elasticity(Mesh*);
//    void setDirichlet(Eigen::VectorXd& v,
//        std::vector<int> dirInd);
//    void setDirichlet(SpMat &mat,
//        std::vector<int> dirInd);
//    SpMat& getK(){return m_K;}
//    Eigen::VectorXd& getRHS(){return m_rhs;}

  private:
    void howMuchAssembled(int,int);

};
