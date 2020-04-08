#pragma once
#include <vector>
#include "types.hpp"


class Domain;

class InterfaceOperatorG
{


  public:
    InterfaceOperatorG(Domain* _domain);
    ~InterfaceOperatorG(){}
//    void multBt(Eigen::MatrixXd& in, Eigen::MatrixXd& out);
//    void multB(Eigen::MatrixXd& in, Eigen::MatrixXd& out);
    void FetiCoarseSpace(std::vector<int>&);
  private:

    Domain* m_p_domain;
    std::vector<int> m_ibuf;
    std::vector<double> m_dbufToSend;
    std::vector<double> m_dbufToRecv;

    std::vector<int> m_listOfNeighbours;
    std::vector<int> m_listOfNeighboursColumPtr;

    void _placeBlockInGlobalGtG(
        std::vector<int>& I_COO,
        std::vector<int>& J_COO,
        std::vector<double>& V_COO,
//        std::vector<T>& _tr,
        Eigen::Ref<Eigen::MatrixXd> _GtG,
        int& tripletOffset, int rowOffset, int colOffset);
    SpMat m_spmatGtG;





};
