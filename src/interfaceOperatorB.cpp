#include "../include/interfaceOperatorB.hpp"
#include "../include/domain.hpp"
#include <iostream>




InterfaceOperatorB::InterfaceOperatorB(Domain* p_dom)
{
  m_p_domain = p_dom;


//  std::cout << "my rank is: " << _dom->GetRank() << std::endl;
//  std::cout << " and my interfaces ranks are: " << std::endl;
//  for (auto& ii : interfaces)
//  {
//    std::cout<<  " " << ii.GetNeighbRank() << std::endl;
//  }
}




void InterfaceOperatorB::multBt(Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  int nRHS = in.cols();
  out.resize(m_p_domain->GetNumberOfPrimalDOFs(),nRHS);
  out.setZero();

  auto interfaces = m_p_domain->GetInterfaces();

  int offset = 0;

  for (auto& iItf : interfaces)
  {
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int nDOFs = (int) iItf.m_interfaceDOFs.size();
    for (int row = 0; row < nDOFs;row++)
      out.row(iItf.m_interfaceDOFs[row]) += scale * in.row(offset + row);
    offset += nDOFs;
  }

}


void InterfaceOperatorB::multB(Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  int nRHS = in.cols();
  out.resize(m_p_domain->GetNumberOfDualDOFs(), nRHS);
  out.setZero();

  auto interfaces = m_p_domain->GetInterfaces();

  int offset = 0;

  for (auto& iItf : interfaces)
  {
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int nDOFs = (int) iItf.m_interfaceDOFs.size();

    m_dbufToRecv.resize(nDOFs * nRHS);
    m_dbufToSend.resize(nDOFs * nRHS);

    for (int row = 0; row < nDOFs; row++)
    {
      out.row(row + offset) = scale * in.row(iItf.m_interfaceDOFs[row]);
      for (int col = 0; col < nRHS; col++)
      {
        m_dbufToSend[nRHS * row + col] = 
          scale * in(iItf.m_interfaceDOFs[row],col);
      }
    }

    m_p_domain->hmpi.SendInt(m_dbufToSend.data(),nDOFs * nRHS, iItf.GetNeighbRank());
    m_p_domain->hmpi.RecvInt(m_dbufToRecv.data(),nDOFs * nRHS, iItf.GetNeighbRank());

    for (int row = 0; row < nDOFs; row++)
    {
      for (int col = 0; col < nRHS; col++)
      {
        out(row + offset, col) += m_dbufToRecv[nRHS * row + col];
      }
    }


    offset += nDOFs;

  }

}
