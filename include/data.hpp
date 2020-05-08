#pragma once

#include <mpi.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "types.hpp"
#include "domain.hpp"
#include "interfaceOperatorB.hpp"
#include "interfaceOperatorG.hpp"
#include "solver.hpp"
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <iostream>
#include <fstream>


class Data{
  public:
    Data(MPI_Comm*);
    ~Data();

    void SymbolicAssembling(std::vector<int>& glbIds);
    void SymbolicAssembling(int[],int size);
    int FinalizeSymbolicAssembling();

    void NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS);
    void NumericAssembling(int glbIds[],
        double valLocK[], double valLocRHS[], int nDofs);

    void FinalizeNumericAssembling();

    Domain* GetDomain(){return  &m_domain;}
    void SetDirichletDOFs(std::vector<int>&v);
    void SetDirichletDOFs(int[],int);
    InterfaceOperatorB* GetInterfaceOperatorB(){return m_p_interfaceOperatorB;}
    void Solve(Eigen::Ref<Eigen::MatrixXd>);
    void Solve(double vals[], int nrows, int ncols = 1);
    void Finalize();
    void PathToSolverOptionFile(std::string path2file);
    boost::property_tree::ptree GetChild(std::string _s0){return m_root.get_child(_s0);}


  private:
    MPI_Comm m_comm;
    int m_mpiRank;
    int m_mpiSize;

    std::vector<int> m_defectPerSubdomains;
    std::vector<int> m_container;
    std::vector<int> m_DirichletGlbDofs;
    Domain m_domain;
    InterfaceOperatorB* m_p_interfaceOperatorB;

    int m_cnt_setLocalMatrix;

    void m_SetKernelNumbering();
    std::vector<int>& _GetDefectPerSubdomains(){return m_defectPerSubdomains;}
    Solver m_solver;


    std::streambuf* m_p_sbuf;
    std::streambuf* m_p_backup;
    std::ofstream m_filestr;
  
    //depracted
    void _dumpTxtFiles(Eigen::Ref<Eigen::MatrixXd>);
    boost::property_tree::ptree m_root;
    int m_verboseLevel;






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
    void m_dbg_printStats();

};
