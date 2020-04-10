#pragma once
#include <iostream>
#include <vector>
#include <map>
#include "types.hpp"
#include <string>
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>



class StiffnessMatrix
{

  public:
    StiffnessMatrix(std::map<int,int>* pg2l,int rank);
    ~StiffnessMatrix(){}


    void AddElementContribution(std::vector<int>& glbIds,
        std::vector<double>& valLocK);
    void AddRHSContribution(std::vector<int>& glbIds,
        std::vector<double>& valLocF);
    void FinalizeNumericPart(const std::vector<int>& DirDOFs);
    void solve(const Eigen::MatrixXd& in, Eigen::MatrixXd& out);
    Eigen::MatrixXd* GetKernel(){ return &m_kerK;}
    int GetDefect(){return m_kerK.cols();}
    void ApplyDirichletBC(SpMat& spmat, std::vector<int> dirInd);
    const Eigen::VectorXd& GetRHS()const {return  m_rhs;};
    void Mult(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void Precond(const Eigen::MatrixXd&, Eigen::MatrixXd&);

  private:
    Eigen::MatrixXd m_kerK;
    std::vector<int> m_nullPivots;
    Eigen::PardisoLU <SpMat> m_pardisoSolver;

    int m_cnt_setLocalMatrix;
    int m_cnt_setLocalRHS;
    int m_numberOfElements;
    std::vector<T> m_trK;
    std::map<int,int>* m_pg2l;
    SpMat m_spmatK;
    int m_myrank;
    Eigen::VectorXd m_rhs;

    int m_neq;

    Eigen::MatrixXd _GetKernelFromK(const SpMat& K);

    void m_FactorizeLinearOperator(std::vector<int> = std::vector<int>(0));
    std::vector<int> GetNullPivots(const Eigen::MatrixXd &);

//DBG ##########
//
  void m_dbg_printStiffnessMatrix(const SpMat&,std::string="");
  void m_dbg_printStiffnessMatrixSingularValues();
  PRECONDITIONER_TYPE m_preconditionerType;
  void _SetDirichletPrecond();

};
