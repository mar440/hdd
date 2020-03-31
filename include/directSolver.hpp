#pragma once

#include "mkl_pardiso.h"
#include "mkl_types.h"
#include <vector>
#include "types.hpp"



using namespace std;
using namespace Eigen;

#include <iostream> 

class DirectSolver{
  public:
    DirectSolver();
    ~DirectSolver();
    void num_factor();
    void sym_factor(SpMat &);
    void solve(Eigen::VectorXd&, Eigen::VectorXd&);

    int testUnsym();
    int testUnsym2();


  private: 
    MKL_INT iparm[64];
    void *pt[64];
    MKL_INT maxfct, mnum, phase, error, msglvl;
    MKL_INT nrhs;
    MKL_INT mtype;
    MKL_INT nrows;
    MKL_INT idum;
    double  ddum;
    MKL_INT nnz;
    vector < double > valuePtr;
    vector < int > innerIndexPtr;
    vector < int > outerIndexPtr;


};







