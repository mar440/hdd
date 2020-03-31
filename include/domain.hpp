#pragma once
#include <vector>
#include <map>
#include <mpi.h>

#include "interface.hpp"
#include "types.hpp"
#include "hmpi.hpp"


class Domain 
{

  public:
//    Domain(int);
    Domain(MPI_Comm*);
    ~Domain(){}

    void SetNeighboursRanks(const std::vector<int>&) ;
    void SetMappingLoc2Glb(std::vector<int>&);

    void NumericAssemblingStiffnessAndRhs(
        std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS);
    void FinalizeStiffnessMatrixAndRhs();
    void SetInterfaces();
    Hmpi hmpi;

  private:

    void m_Init();
    int m_rank;
    MPI_Comm* m_pcomm;

    void m_SetMappingGlb2Loc();
    std::vector<int> m_l2g;
    std::map<int,int> m_g2l;

    std::vector<int> m_neighboursRanks;

    int m_neq;
    int m_l0;
    int m_cnt_setLocalMatrix;

    // A * x = b
    std::vector<T> m_trK;
    SpMat m_stiffnessMatrix;
    Eigen::VectorXd m_rhs;

    //
    std::vector<Interface> m_interfaces;





    // dbg - 
    void m_dbg_printNeighboursRanks();
    void m_dbg_print_l2g();
    void m_dbg_printStiffnessMatrix();
    void m_dbg_printStiffnessMatrixSingularValues();








};
