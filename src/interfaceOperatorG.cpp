#if 0
#include "../include/interfaceOperatorG.hpp"

#include "../include/domain.hpp"
#include "../include/stiffnessMatrix.hpp"
#include "../include/linearAlgebra.hpp"
#include <iostream>
#include <Eigen/Dense>
#include <cmath>




InterfaceOperatorG::InterfaceOperatorG(Domain* p_dom)
{
  m_p_domain = p_dom;

  m_listOfNeighbours.resize(0);
  m_listOfNeighboursColumPtr.resize(0);
  m_spmatGtG.resize(0,0);
  m_p_defectPerSubdomains = nullptr;

}



void InterfaceOperatorG::_placeBlockInGlobalGtG(
    std::vector<int>& I_COO,
    std::vector<int>& J_COO,
    std::vector<double>& V_COO,
    Eigen::Ref<Eigen::MatrixXd> _GtG,
    int& tripletOffset, int offsetRow, int offsetCol)
{
    for (int rowG = 0; rowG <  _GtG.rows();rowG++)
    {
      for (int colG = 0; colG <  _GtG.cols();colG++)
      {
        int _row = rowG + offsetRow;
        int _col = colG + offsetCol;
        double _val = _GtG(rowG,colG);
//        _tr[tripletOffset] = T(_row,_col,_val);

        I_COO[tripletOffset] = _row;
        J_COO[tripletOffset] = _col;
        V_COO[tripletOffset] = _val;

        tripletOffset++;
      }
    }
}



void InterfaceOperatorG::FetiCoarseSpace(
    std::vector<int>& defectPerSubdomains)
{

  m_p_defectPerSubdomains =  &defectPerSubdomains;
  _FetiCoarseSpaceAssembling();

  m_pardisoSolver.analyzePattern(m_spmatGtG);
  m_pardisoSolver.factorize(m_spmatGtG);

}

void InterfaceOperatorG::solve(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  out = m_pardisoSolver.solve(in);

}




void InterfaceOperatorG::_FetiCoarseSpaceAssembling()
{

  auto kerK = m_p_domain->GetStiffnessMatrix()->GetKernel();
  int myDefect = m_p_domain->GetStiffnessMatrix()->GetDefect(); // = kerK.cols()
  auto interfaces = m_p_domain->GetInterfaces();

  // buffers
  std::vector<Eigen::MatrixXd> BR_myRank(interfaces.size(),Eigen::MatrixXd(0,0));
  std::vector<Eigen::MatrixXd> BR_neighb(interfaces.size(),Eigen::MatrixXd(0,0));

  // GtG
  //  - diagonal block
  Eigen::MatrixXd GtG_local_diag = Eigen::MatrixXd::Zero(myDefect,myDefect);
  //  - offdiagonal blocks 
  std::vector<Eigen::MatrixXd> GtG_local_offdiag(interfaces.size(),Eigen::MatrixXd(0,0));


  int nRHS_myRank = myDefect;


  // each rank computes whole block row 
  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++)
  {

    Interface &iItf = interfaces[cntI];

    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int neqIntf = (int) iItf.m_interfaceDOFs.size();

    int nRHS_neighb = (*m_p_defectPerSubdomains)[iItf.GetNeighbRank()];


    BR_myRank[cntI].resize(neqIntf,nRHS_myRank);
    BR_neighb[cntI].resize(neqIntf,nRHS_neighb);

    for (int row = 0; row < neqIntf; row++)
      BR_myRank[cntI].row(row) = scale * (*kerK).row(iItf.m_interfaceDOFs[row]);

    m_p_domain->hmpi.SendDbl(BR_myRank[cntI].data(),neqIntf * nRHS_myRank, iItf.GetNeighbRank());
    m_p_domain->hmpi.RecvDbl(BR_neighb[cntI].data(),neqIntf * nRHS_neighb, iItf.GetNeighbRank());

#if DBG > 4
    std::cout<< "BR_myRank - myRank: " << m_p_domain->GetRank() << std::endl;
    std::cout<< BR_myRank[cntI]   << std::endl;

    std::cout<< "BR_neighb - neighb: " << iItf.GetNeighbRank() << std::endl;
    std::cout<< BR_neighb[cntI]   << std::endl;
#endif
    // diagonal GtG(i = myRank, j = myRank)
    GtG_local_diag += BR_myRank[cntI].transpose() * BR_myRank[cntI];
    // offdiagonal GtG(i = myRank, j = neighbours)
    GtG_local_offdiag[cntI] = BR_myRank[cntI].transpose() * BR_neighb[cntI];
  }

  // empty buffers
  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++){
    BR_myRank[cntI].resize(0,0); BR_neighb[cntI].resize(0,0);
  }


#if DBG > 3
  std::cout<< "GtG_local_diag: \n";

  std::cout <<  '[' << m_p_domain->GetRank() <<',' << m_p_domain->GetRank()  << ']' << '\n';
  std::cout<< GtG_local_diag  << std::endl;

  for (int cntI = 0; cntI < (int) GtG_local_offdiag.size(); cntI++)
  {
    std::cout <<  '[' << m_p_domain->GetRank() <<',' <<interfaces[cntI].GetNeighbRank()  << ']' << '\n';
    std::cout<< GtG_local_offdiag[cntI]   << std::endl;
  }
#endif
  
  // GtG distributed per ranks 


  // each rank counts number of interfaces
  int numberOfInterfacesPerRank = (int)interfaces.size();

  //
  const int root = 0;
  std::vector<int> numberOfNeighboursRoot(0);
  if (m_p_domain->GetRank() == 0)
      numberOfNeighboursRoot.resize(m_p_domain->GetNumberOfSubdomains(),-1);
  std::cout << "numberOfNeighboursRoot: << "  << numberOfNeighboursRoot.size() << '\n';

  m_p_domain->hmpi.GatherInt(&numberOfInterfacesPerRank, 1,
               numberOfNeighboursRoot.data(), 1 ,root);


  if (m_p_domain->GetRank() == 0)
  for (auto& ii : numberOfNeighboursRoot) std::cout << ii << ' ';
  std::cout << '\n';

  int sumOfInterfacesGlobal(0);
  m_p_domain->hmpi.ReduceInt(&numberOfInterfacesPerRank,&sumOfInterfacesGlobal,
      1, MPI_SUM, root);

  std::cout <<  "numberOfInterfacesPerRank: " << numberOfInterfacesPerRank<< std::endl;
  if (m_p_domain->GetRank()  == root)
    std::cout <<  "sumOfInterfacesGlobal: " << sumOfInterfacesGlobal << std::endl;



  std::vector<int> listOfInterfacesLocal(0);
  for (auto& itf : interfaces)
    listOfInterfacesLocal.push_back(itf.GetNeighbRank());


  int edgeSum=0;
  if (m_p_domain->GetRank()  == root){
    m_listOfNeighbours.resize(sumOfInterfacesGlobal,-1);

    m_listOfNeighboursColumPtr.resize(m_p_domain->GetNumberOfSubdomains()+1);
    for(int i=0; i<m_p_domain->GetNumberOfSubdomains(); ++i) {
        m_listOfNeighboursColumPtr[i]=edgeSum;
        edgeSum+=numberOfNeighboursRoot[i];
    }
    m_listOfNeighboursColumPtr[m_p_domain->GetNumberOfSubdomains()]=edgeSum;
  }



    m_p_domain->hmpi.GathervInt(
        listOfInterfacesLocal.data(), listOfInterfacesLocal.size(), // to send
        m_listOfNeighbours.data(), numberOfNeighboursRoot.data(),   // info
        m_listOfNeighboursColumPtr.data(),  root);


  if (m_p_domain->GetRank()  == root){
      std::cout <<  "===========================+\n";


      for (int in = 0; in < (int) numberOfNeighboursRoot.size(); in++ )
      {

        std::cout << "(" << in << "): " ;
        for (int jn = m_listOfNeighboursColumPtr[in];
            jn < m_listOfNeighboursColumPtr[in+1]; jn++)
        {
          std::cout << m_listOfNeighbours[jn]<< ' ';
        }
        std::cout << '\n';
      }
      std::cout <<  "===========================+\n";


  }


  //
  // send blocks to root
  //


  int sizeOfSendingBuffer = myDefect * myDefect;

  int cntI(0);
  for (auto& itf : interfaces){
    int tmp0 = (*m_p_defectPerSubdomains)[itf.GetNeighbRank()] * myDefect;
    int tmp1 = GtG_local_offdiag[cntI++].size();
    sizeOfSendingBuffer += tmp0;
    std::cout << "tmp0 & tmp1: " << tmp0 << ' ' << tmp1 << '\n';
  }

  std::cout << "sizeOfSendingBuffer: " << sizeOfSendingBuffer<< '\n';


  std::vector<double> sendingBuffer(sizeOfSendingBuffer,0);

  // put diagonal block in buffer
  
//  sendingBuffer
  int offset(0);
  std::memcpy( sendingBuffer.data(),GtG_local_diag.data(),
      sizeof( double ) * GtG_local_diag.size());

  offset += GtG_local_diag.size();

  for (auto& GtG_ij : GtG_local_offdiag)
  {
    std::memcpy( sendingBuffer.data() + offset,GtG_ij.data(),
      sizeof( double ) * GtG_ij.size());
    offset += GtG_ij.size();
  }

#if DBG > 3
  for (auto& ii : sendingBuffer){std::cout << ii << ' ';} std::cout << '\n';
#endif


  std::vector<int> sizePerRank(0);
  std::vector<double> receivingBuffer(0);
  std::vector<int> sizePerRankColumPtr(0);

  if (m_p_domain->GetRank()  == root){



    int cntAllElements(0);
    for (int in = 0; in < (int) numberOfNeighboursRoot.size(); in++ )
    {
      int diagRank = in;

      int defectDiag = (*m_p_defectPerSubdomains)[diagRank];

      sizePerRank.push_back(pow(defectDiag,2));

      sizePerRankColumPtr.push_back(cntAllElements);
      for (int jn = m_listOfNeighboursColumPtr[in];
          jn < m_listOfNeighboursColumPtr[in+1]; jn++)
      {

        int offdiagRank = m_listOfNeighbours[jn];
        int defectOffdiag = (*m_p_defectPerSubdomains)[offdiagRank];

        sizePerRank.back() += defectDiag * defectOffdiag;
      }
      cntAllElements += sizePerRank.back();
    }
    sizePerRankColumPtr.push_back(cntAllElements);


    receivingBuffer.resize(cntAllElements,-1);
  }

  std::cout << "receivingBuffer.size() = " << receivingBuffer.size() << '\n';



  std::cout << "sendingBuffer.size(): " << sendingBuffer.size() << '\n';

  m_p_domain->hmpi.GathervDbl(
    sendingBuffer.data(), sendingBuffer.size(),
    receivingBuffer.data(), sizePerRank.data(),
    sizePerRankColumPtr.data(),  root);

//  for (int ii = 0; ii < (int) sizePerRank.size(); ii++)
//  {
//    for (int jj = sizePerRankColumPtr[ii];jj < sizePerRankColumPtr[ii+1];jj++) 
//    {
//      std::cout << receivingBuffer[jj] << ' ';
//    }
//    std::cout << '\n';
//  }

  std::vector<int>    I_COO_GtG(0);
  std::vector<int>    J_COO_GtG(0);
  std::vector<double> V_COO_GtG(0);


  int defectGtG(0);
  int nnzGtG(0);

  if (m_p_domain->GetRank()  == root)
  {
    nnzGtG = (int)receivingBuffer.size();

    for (auto& ii : (*m_p_defectPerSubdomains))
      defectGtG += ii;
  }

  int tmpInt[] = {nnzGtG, defectGtG};

  m_p_domain->hmpi.BcastInt(&tmpInt,2,root);

  nnzGtG = tmpInt[0];
  defectGtG = tmpInt[1];


  std::cout << "defectGtG:  " << defectGtG << '\n'; 
  std::cout << "nnzGtG:     " << nnzGtG << '\n'; 


  if (m_p_domain->GetRank()  == root)
  {

    I_COO_GtG.resize(nnzGtG);
    J_COO_GtG.resize(nnzGtG);
    V_COO_GtG.resize(nnzGtG);


    std::vector<int>  cumulativeDefectPerSubdomains((*m_p_defectPerSubdomains).size(),0);
    for (int ii = 0; ii < (int) (*m_p_defectPerSubdomains).size() - 1 ; ii++)
    {
      cumulativeDefectPerSubdomains[ii+1] = 
        cumulativeDefectPerSubdomains[ii] + (*m_p_defectPerSubdomains)[ii];
    }

    int cntGtG(0);
    int offsetRow(0);
    int offsetCol(0);
    int cntAllElements(0);
    for (int in = 0; in < (int) numberOfNeighboursRoot.size(); in++ )
    {
      int diagRank = in;
      int defectDiag = (*m_p_defectPerSubdomains)[diagRank];

      Eigen::Map<Eigen::MatrixXd> 
        GtG_local_diag(receivingBuffer.data() + cntAllElements,
          defectDiag, defectDiag);
      cntAllElements +=  defectDiag * defectDiag;

      offsetRow = cumulativeDefectPerSubdomains[diagRank];
      offsetCol = cumulativeDefectPerSubdomains[diagRank];

      _placeBlockInGlobalGtG(I_COO_GtG,J_COO_GtG,V_COO_GtG,
       //   trGtG, 
       GtG_local_diag,cntGtG,offsetRow ,offsetCol);



#if DBG > 3
      std::cout <<  '[' << in <<',' << in << ']' << '\n';
      std::cout << GtG_local_diag << '\n';
      std::cout << '\n';
#endif

      for (int jn = m_listOfNeighboursColumPtr[in];
          jn < m_listOfNeighboursColumPtr[in+1]; jn++)
      {

        int offdiagRank = m_listOfNeighbours[jn];
        int defectOffdiag = (*m_p_defectPerSubdomains)[offdiagRank];

        Eigen::Map<Eigen::MatrixXd> 
          GtG_local_offdiag(receivingBuffer.data() + cntAllElements,
            defectDiag, defectOffdiag);


        offsetRow = cumulativeDefectPerSubdomains[diagRank];
        offsetCol = cumulativeDefectPerSubdomains[offdiagRank];

        _placeBlockInGlobalGtG(I_COO_GtG,J_COO_GtG,V_COO_GtG,
           // trGtG, 
            GtG_local_offdiag,cntGtG,offsetRow, offsetCol);

#if DBG > 2
        std::cout <<  '[' << in <<',' << offdiagRank << ']' << '\n';
        std::cout << GtG_local_offdiag << '\n';
        std::cout << '\n';
#endif

        cntAllElements += defectDiag * defectOffdiag;
      }
    }
  }

  I_COO_GtG.resize(nnzGtG);
  J_COO_GtG.resize(nnzGtG);
  V_COO_GtG.resize(nnzGtG);
//
  m_p_domain->hmpi.BcastInt(I_COO_GtG.data(),nnzGtG,root);
  m_p_domain->hmpi.BcastInt(J_COO_GtG.data(),nnzGtG,root);
  m_p_domain->hmpi.BcastDbl(V_COO_GtG.data(),nnzGtG,root);

  std::vector<T> trGtG(nnzGtG,T(0,0,0));

  for (int iG = 0; iG < nnzGtG; iG++)
    trGtG[iG] = T(I_COO_GtG[iG],J_COO_GtG[iG],V_COO_GtG[iG]);

  m_spmatGtG.resize(defectGtG, defectGtG);
  m_spmatGtG.setFromTriplets(trGtG.begin(),trGtG.end());
//
  std::string fname = "GtG"; fname += std::to_string(m_p_domain->GetRank()) + ".txt";
  tools::printMatrix(m_spmatGtG,fname);


}
#endif
