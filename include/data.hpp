#pragma once

#include <mpi.h>
#include <Eigen/Sparse>
#include <Eigen/Dense>

#include "mesh.hpp"
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
    void FinalizeSymbolicAssembling();

    void NumericAssembling(std::vector<int>& glbIds,
        std::vector<double>& valLocK,
        std::vector<double>& valLocRHS);
    void FinalizeNumericAssembling();

    Domain* GetDomain(){return  &m_domain;}
    void SetDirichletDOFs(std::vector<int>&v);
    InterfaceOperatorB* GetInterfaceOperatorB(){return m_p_interfaceOperatorB;}
    void Solve(Eigen::VectorXd&);
    void Finalize();
    void ParseJsonFile(std::string path2file);
//    InterfaceOperatorG* GetInterfaceOperatorG(){return m_p_interfaceOperatorG;}
    boost::property_tree::ptree GetChild(std::string _s0){return m_root.get_child(_s0);}
//


  private:
    MPI_Comm* m_pcomm;
    int m_mpiRank;
    int m_mpiSize;

    std::vector<int> m_defectPerSubdomains;
    std::vector<int> m_container;
    std::vector<int> m_DirichletGlbDofs;
    Domain m_domain;
    InterfaceOperatorB* m_p_interfaceOperatorB;
    //InterfaceOperatorG* m_p_interfaceOperatorG;

    int m_cnt_setLocalMatrix;

    void m_SetKernelNumbering();
    std::vector<int>& _GetDefectPerSubdomains(){return m_defectPerSubdomains;}
    Solver m_solver;


    std::streambuf* m_p_sbuf;
    std::streambuf* m_p_backup;
    std::ofstream m_filestr;

    void _dumpTxtFiles(Eigen::VectorXd&);
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
