#pragma once
#include <vector>
#include <Eigen/Dense>


class Domain;

class InterfaceOperatorB
{


  public:
    InterfaceOperatorB(Domain* _domain);
    ~InterfaceOperatorB(){}
    void multBt(Eigen::MatrixXd& in, Eigen::MatrixXd& out);
    void multB(Eigen::MatrixXd& in, Eigen::MatrixXd& out);
  private:

    Domain* m_p_domain;
    std::vector<int> m_ibuf;
    std::vector<double> m_dbufToSend;
    std::vector<double> m_dbufToRecv;






};
