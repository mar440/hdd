#pragma once
#include <vector>
#include <Eigen/Dense>
#include "types.hpp"
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>


class Domain;

class InterfaceOperatorB
{


  public:
    InterfaceOperatorB(Domain* _domain);
    ~InterfaceOperatorB(){}
    void multBt(const Eigen::MatrixXd& in, Eigen::MatrixXd& out);
    void multB(const Eigen::MatrixXd& in, Eigen::MatrixXd& out);
  // orig G
    void FetiCoarseSpace(std::vector<int>&);

//    void GinvGtG(const Eigen::MatrixXd& in, Eigen::MatrixXd& lambda0);

    void mult_invGtG(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    Eigen::MatrixXd Projection(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void mult_F(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void Scaling(const Eigen::MatrixXd&, Eigen::MatrixXd&);
    void Scaling(Eigen::MatrixXd&);
  private:

    Domain* m_p_domain;
    int m_root;
    Eigen::MatrixXd m_dbufToSend;
    Eigen::MatrixXd m_dbufToRecv;

    // orig G
    void _FetiCoarseSpaceAssembling();
    std::vector<int>* m_p_defectPerSubdomains;


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
    int m_GtG_dim;

    std::vector<int> m_cumulativeDefectPerSubdomains;
    std::vector<int> m_numberOfNeighboursRoot;


    Eigen::PardisoLU<SpMat> m_pardisoSolver;

    void _SetScaling();

    std::vector<double> m_scaling;

};
