#pragma once
#include <vector>
#include <map>
#include <mpi.h>

#include "interface.hpp"
#include "types.hpp"
#include "hmpi.hpp"



class StiffnessMatrix;


class Domain 
{

  public:
    Domain(MPI_Comm*);
    ~Domain();

    void SetNeighboursRanks(const std::vector<int>&) ;
    void SetMappingLoc2Glb(std::vector<int>&);

    void SetInterfaces();
    void _SetInterfaces();
    const std::vector<Interface>& GetInterfaces() {return m_interfaces;};
    Hmpi hmpi;
    int GetRank(){return m_mpirank;}
    int GetMpiSize(){return m_mpisize;}
    int GetNumberOfSubdomains(){return m_mpisize;}
    int GetNumberOfPrimalDOFs(){return m_neqPrimal;}
    int GetNumberOfDualDOFs(){return m_neqDual;}
    void InitStiffnessMatrix(PRECONDITIONER_TYPE);
    StiffnessMatrix* GetStiffnessMatrix(){return m_p_stiffnessMatrix;}
    void SetDirichletDOFs(std::vector<int>& glbDirDOFs);
    const std::vector<int>&  GetDirichletDOFs()const{return m_DirichletDOFs;}
    std::vector<double>* GetMultiplicity(){return &m_multiplicity;}
    void HandlePreconditioning();


  private:

    void m_Init();
    int m_mpirank;
    int m_mpisize;
    MPI_Comm m_comm;

    void m_SetMappingGlb2Loc();
    std::vector<int> m_l2g;
    std::map<int,int> m_g2l;

    std::vector<int> m_neighboursRanks;

    int m_neqPrimal;
    int m_neqDual;
    int m_l0;
    std::vector<int> m_DirichletDOFs;

    std::vector<int> m_I_DirichletPrecondDOFs;
    std::vector<int> m_B_DirichletPrecondDOFs;

    // A * x = b
    StiffnessMatrix* m_p_stiffnessMatrix;
    Eigen::VectorXd m_rhs;
    Eigen::MatrixXd m_matG;

    // for interface operator
    std::vector<Interface> m_interfaces;
    std::map<int,int> m_intf_g2l;

    std::vector<double> m_multiplicity;

    void _SetDirichletPrecondDOFs();
    void _SearchNeighbours(std::vector<int>&);

    int m_verboseLevel;

    // dbg ---------------------------- 
    void m_dbg_printNeighboursRanks();
    void m_dbg_print_l2g();

    std::vector<int> m_listOfNeighbours;
    std::vector<int> m_listOfNeighboursColumPtr;
};
