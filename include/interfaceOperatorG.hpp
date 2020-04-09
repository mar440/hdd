#if 0
#pragma once
#include <vector>
#include "types.hpp"
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>

class Domain;
class InterfaceOperatorG
{

  public:
    InterfaceOperatorG(Domain* _domain);
    ~InterfaceOperatorG(){}
    void FetiCoarseSpace(std::vector<int>&);
  private:

    void _FetiCoarseSpaceAssembling();
    std::vector<int>* m_p_defectPerSubdomains;


    Domain* m_p_domain;
    std::vector<double> m_dbufToSend;
    std::vector<double> m_dbufToRecv;

    std::vector<int> m_listOfNeighbours;
    std::vector<int> m_listOfNeighboursColumPtr;
    void solve(const Eigen::MatrixXd&, Eigen::MatrixXd&);

    void _placeBlockInGlobalGtG(
        std::vector<int>& I_COO,
        std::vector<int>& J_COO,
        std::vector<double>& V_COO,
//        std::vector<T>& _tr,
        Eigen::Ref<Eigen::MatrixXd> _GtG,
        int& tripletOffset, int rowOffset, int colOffset);
    SpMat m_spmatGtG;

    Eigen::PardisoLU<SpMat> m_pardisoSolver;





};
#endif
