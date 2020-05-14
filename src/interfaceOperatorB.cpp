#include "../include/interfaceOperatorB.hpp"
#include "../include/domain.hpp"
#include <iostream>


#include "../include/stiffnessMatrix.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/linearAlgebra.hpp"
#include <cmath>



#define GTG_ALL_NODES
#define TAG0 401
#define TAG1 402

InterfaceOperatorB::InterfaceOperatorB(Domain* p_dom)
{

  m_p_domain = p_dom;
  m_listOfNeighbours.resize(0);
  m_listOfNeighboursColumPtr.resize(0);
  m_spmatGtG.resize(0,0);
  m_p_defectPerSubdomains = nullptr;
  m_root = 0;
  m_cumulativeDefectPerSubdomains.resize(0);
  m_numberOfNeighboursRoot.resize(0);
  m_GtG_dim = 0;
  m_dbufToSend.resize(0,0);
  m_dbufToRecv.resize(0,0);
  m_scaling.resize(0);
  _SetScaling();

}




void InterfaceOperatorB::multBt(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  int nRHS = in.cols();
  out.resize(m_p_domain->GetNumberOfPrimalDOFs(),nRHS);
  out.setZero();

  auto interfaces = m_p_domain->GetInterfaces();

//  int offset = 0;

  for (auto& iItf : interfaces)
  {
    ////////////////////////////////////////////////////
    // CONVENTION: subdomain with higher rank has
    // operator B (or G) with negative coefficients (-1)
    ////////////////////////////////////////////////////
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int nDOFs = (int) iItf.m_interfaceDOFs.size();
    for (int row = 0; row < nDOFs;row++)
      out.row(iItf.m_interfaceDOFs[row]) += scale * in.row(iItf.m_offset + row);
  }
}

void InterfaceOperatorB::multB(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{


  int nRHS = in.cols();
  out.resize(m_p_domain->GetNumberOfDualDOFs(), nRHS);
  out.setZero();

  auto interfaces = m_p_domain->GetInterfaces();


  int nInterf = (int)interfaces.size();
  MPI_Request requests[nInterf];
  MPI_Status statuses[nInterf];


  m_sendBuffers.resize(nInterf,Eigen::MatrixXd::Zero(0,0));
  for (int intfId = 0; intfId < nInterf; intfId++)
  {
    Interface& iItf = interfaces[intfId];
    ////////////////////////////////////////////////////
    // CONVENTION: subdomain with higher rank has
    // operator B (or G) with negative coefficients (-1)
    ////////////////////////////////////////////////////
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int nDOFs = (int) iItf.m_interfaceDOFs.size();

    m_sendBuffers[intfId].resize(nDOFs , nRHS);

    for (int row = 0; row < nDOFs; row++)
    {
      out.row(row + iItf.m_offset) += scale * in.row(iItf.m_interfaceDOFs[row]);
      m_sendBuffers[intfId].row(row) = scale * in.row(iItf.m_interfaceDOFs[row]);
    }
    MPI_Isend(m_sendBuffers[intfId].data(),nDOFs * nRHS, MPI_DOUBLE, iItf.GetNeighbRank(), TAG0,
      m_p_domain->hmpi.GetComm(), &requests[intfId]);
  }

  HDDTRACES
  int msg_avail = -1;
  MPI_Status status1, status2;
  int bufersize(0);
  int cntRecvMsg(0);

  while(true)
  {

    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_p_domain->hmpi.GetComm(),&msg_avail, &status1);

    HDDTRACES
    if (msg_avail!=0)
    {
      HDDTRACES
      MPI_Get_count(&status1, MPI_DOUBLE, &bufersize);

      int neighRank = status1.MPI_SOURCE;

      auto intf_g2l = m_p_domain->GetInterfacesMapping_g2l();
      int intfId = intf_g2l[neighRank];

      Interface& iItf = interfaces[intfId];

      int nDOFs = (int) iItf.m_interfaceDOFs.size();
      m_dbufToRecv.resize(nDOFs , nRHS);

//      std::cout << "same numbers: " << bufersize << " " << nDOFs * nRHS << '\n';

      MPI_Recv(m_dbufToRecv.data(),nDOFs * nRHS, MPI_DOUBLE,  iItf.GetNeighbRank(),TAG0,
        m_p_domain->hmpi.GetComm(), &status2);


      for (int row = 0; row < nDOFs; row++)
        out.row(row + iItf.m_offset) += m_dbufToRecv.row(row);

      if (++cntRecvMsg == nInterf) break;
    }
  }
  MPI_Wait(requests,statuses);
  m_p_domain->hmpi.Barrier();

  m_releaseBuffers();


}


void InterfaceOperatorB::_solve(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

#if !defined(GTG_ALL_NODES)
  if (m_p_domain->GetRank()  == m_root)
#endif
    out = m_pardisoSolver.solve(in);

}


void InterfaceOperatorB::_placeBlockInGlobalGtG(
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

      I_COO[tripletOffset] = _row;
      J_COO[tripletOffset] = _col;
      V_COO[tripletOffset] = _val;

      tripletOffset++;
    }
  }
}

void InterfaceOperatorB::FetiCoarseSpace(
    std::vector<int>& defectPerSubdomains)
{

  HDDTRACES
  m_p_defectPerSubdomains =  &defectPerSubdomains;
  HDDTRACES
  _FetiCoarseSpaceAssembling();
  HDDTRACES
#if !defined(GTG_ALL_NODES)
  if (m_p_domain->GetRank()  == m_root) {
#endif
    m_pardisoSolver.analyzePattern(m_spmatGtG);
    m_pardisoSolver.factorize(m_spmatGtG);
#if !defined(GTG_ALL_NODES)
  }
#endif

}


void InterfaceOperatorB::_FetiCoarseSpaceAssembling()
{

  HDDTRACES
  auto kerK = m_p_domain->GetStiffnessMatrix()->GetKernel();
  int myDefect = m_p_domain->GetStiffnessMatrix()->GetDefect(); 
  auto interfaces = m_p_domain->GetInterfaces();

  // buffers
  std::vector<Eigen::MatrixXd> BR_myRank(interfaces.size(),Eigen::MatrixXd(0,0));
  std::vector<Eigen::MatrixXd> BR_neighb(interfaces.size(),Eigen::MatrixXd(0,0));

  HDDTRACES
  // GtG
  //  - diagonal block
  Eigen::MatrixXd GtG_local_diag = Eigen::MatrixXd::Zero(myDefect,myDefect);
  //  - offdiagonal blocks 
  std::vector<Eigen::MatrixXd> GtG_local_offdiag(interfaces.size(),Eigen::MatrixXd(0,0));

  HDDTRACES
  int nRHS_myRank = myDefect;

  int nInterf = (int)interfaces.size();

  MPI_Request requests[nInterf];
  MPI_Status statuses[nInterf];

  HDDTRACES
  // each rank computes whole block row 
  for (int cntI = 0; cntI < nInterf; cntI++)
  {

    Interface &iItf = interfaces[cntI];

    ////////////////////////////////////////////////////
    // CONVENTION: subdomain with higher rank has
    // operator B (or G) with negative coefficients (-1)
    ////////////////////////////////////////////////////
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int neqIntf = (int) iItf.m_interfaceDOFs.size();

    int nRHS_neighb = (*m_p_defectPerSubdomains)[iItf.GetNeighbRank()];


    BR_myRank[cntI].resize(neqIntf,nRHS_myRank);
    BR_neighb[cntI].resize(neqIntf,nRHS_neighb);

    // !!! ALERT in ''Gt*lambda = e'' negative sign: Gt = -Rt*Bt, e = -Rt*f
    for (int row = 0; row < neqIntf; row++)
      BR_myRank[cntI].row(row) = (-1) * scale * (*kerK).row(iItf.m_interfaceDOFs[row]);

// m_p_domain->hmpi.IsendDbl(BR_myRank[cntI].data(),neqIntf * nRHS_myRank, iItf.GetNeighbRank());
    MPI_Isend(BR_myRank[cntI].data(),neqIntf*nRHS_myRank, MPI_DOUBLE, iItf.GetNeighbRank(),TAG0,
      m_p_domain->hmpi.GetComm(), requests);
  }

  int msg_avail = -1;
  MPI_Status status1, status2;
  int bufersize(0);
  int cntRecvMsg(0);

  while(true)
  {

    MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, m_p_domain->hmpi.GetComm(),&msg_avail, &status1);

    HDDTRACES
    if (msg_avail!=0)
    {
      HDDTRACES
      MPI_Get_count(&status1, MPI_DOUBLE, &bufersize);

      int neighRank = status1.MPI_SOURCE;

      auto intf_g2l = m_p_domain->GetInterfacesMapping_g2l();
      int intfId = intf_g2l[neighRank];

      Interface& iItf = interfaces[intfId];
      int nRHS_neighb = (*m_p_defectPerSubdomains)[iItf.GetNeighbRank()];

      int nDOFs = (int) iItf.m_interfaceDOFs.size();
      m_dbufToRecv.resize(nDOFs , nRHS_neighb);

      MPI_Recv(BR_neighb[intfId].data(),nDOFs*nRHS_neighb, MPI_DOUBLE,iItf.GetNeighbRank(),TAG0,
        m_p_domain->hmpi.GetComm(), &status2);

      if (++cntRecvMsg == nInterf) break;
    }
  }


  MPI_Wait(requests,statuses);
  m_p_domain->hmpi.Barrier();

  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++)
  {

//    Interface &iItf = interfaces[cntI];
//
//    int neqIntf = (int) iItf.m_interfaceDOFs.size();
//
//    int nRHS_neighb = (*m_p_defectPerSubdomains)[iItf.GetNeighbRank()];
//
//    m_p_domain->hmpi.IrecvDbl(BR_neighb[cntI].data(),neqIntf * nRHS_neighb, iItf.GetNeighbRank());

#if DBG > 4
    std::cout<< "BR_myRank - myRank: " << m_p_domain->GetRank() << std::endl;
    std::cout<< BR_myRank[cntI]   << std::endl;

    std::cout<< "BR_neighb - neighb: " << iItf.GetNeighbRank() << std::endl;
    std::cout<< BR_neighb[cntI]   << std::endl;
#endif

  }





//  m_p_domain->hmpi.Wait();

  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++)
  {
    // diagonal GtG(i = myRank, j = myRank)
    
    if (BR_myRank[cntI].size() > 0 && BR_myRank[cntI].size() > 0)
    {
      GtG_local_diag += BR_myRank[cntI].transpose() * BR_myRank[cntI];
      // offdiagonal GtG(i = myRank, j = neighbours)
      GtG_local_offdiag[cntI] = BR_myRank[cntI].transpose() * BR_neighb[cntI];
    }
    else
    {
#if DBG > 4
      std::cout<< "no contribution to GtG\n";
#endif
    }
  }
  HDDTRACES

//  m_p_domain->hmpi.Barrier();
  // empty buffers
  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++){
    BR_myRank[cntI].resize(0,0); BR_neighb[cntI].resize(0,0);
  }

  HDDTRACES

#if DBG > 4
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
  std::vector<int> numberOfNeighboursRoot(0);
  if (m_p_domain->GetRank() == 0)
      numberOfNeighboursRoot.resize(m_p_domain->GetNumberOfSubdomains(),-1);
#if DBG > 2
  std::cout << "numberOfNeighboursRoot: << "  << numberOfNeighboursRoot.size() << std::endl;
  std::cout << "numberOfInterfacesPerRank:" << numberOfInterfacesPerRank << std::endl;
#endif

//  m_p_domain->hmpi.Barrier();

  HDDTRACES
  m_p_domain->hmpi.GatherInt(&numberOfInterfacesPerRank, 1,
               numberOfNeighboursRoot.data(), 1 ,m_root);
  HDDTRACES

#if DBG > 2
  if (m_p_domain->GetRank() == 0)
    for (auto& ii : numberOfNeighboursRoot) std::cout << ii << ' ';
  std::cout << '\n';
#endif

  int sumOfInterfacesGlobal(0);
  m_p_domain->hmpi.ReduceInt(&numberOfInterfacesPerRank,&sumOfInterfacesGlobal,
      1, MPI_SUM, m_root);

#if DBG > 2
  std::cout <<  "numberOfInterfacesPerRank: " << numberOfInterfacesPerRank<< std::endl;
  if (m_p_domain->GetRank()  == m_root)
    std::cout <<  "sumOfInterfacesGlobal: " << sumOfInterfacesGlobal << std::endl;
#endif


  HDDTRACES

  std::vector<int> listOfInterfacesLocal(0);
  for (auto& itf : interfaces)
    listOfInterfacesLocal.push_back(itf.GetNeighbRank());

  int edgeSum=0;
  if (m_p_domain->GetRank()  == m_root){
    m_listOfNeighbours.resize(sumOfInterfacesGlobal,-1);

    m_listOfNeighboursColumPtr.resize(m_p_domain->GetNumberOfSubdomains()+1);
    for(int i=0; i<m_p_domain->GetNumberOfSubdomains(); ++i) {
        m_listOfNeighboursColumPtr[i]=edgeSum;
        edgeSum+=numberOfNeighboursRoot[i];
    }
    m_listOfNeighboursColumPtr[m_p_domain->GetNumberOfSubdomains()]=edgeSum;
  }


  HDDTRACES

  m_p_domain->hmpi.GathervInt(
      listOfInterfacesLocal.data(), listOfInterfacesLocal.size(), // to send
      m_listOfNeighbours.data(), 
      numberOfNeighboursRoot.data(),m_listOfNeighboursColumPtr.data(),  m_root);


  HDDTRACES


  int sizeOfSendingBuffer = myDefect * myDefect;

  for (auto& itf : interfaces)
  {
    int tmp0 = (*m_p_defectPerSubdomains)[itf.GetNeighbRank()] * myDefect;
    sizeOfSendingBuffer += tmp0;
  }

  HDDTRACES

#if DBG > 2
  std::cout << "sizeOfSendingBuffer: " << sizeOfSendingBuffer<< '\n';
#endif
  HDDTRACES

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

  if (m_p_domain->GetRank()  == m_root){



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

#if DBG > 3
  std::cout << "receivingBuffer.size() = " << receivingBuffer.size() << '\n';
  std::cout << "sendingBuffer.size(): " << sendingBuffer.size() << '\n';
#endif

  m_p_domain->hmpi.GathervDbl(
    sendingBuffer.data(), sendingBuffer.size(),
    receivingBuffer.data(), sizePerRank.data(),
    sizePerRankColumPtr.data(),  m_root);




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


  int nnzGtG(0);

  if (m_p_domain->GetRank()  == m_root)
  {
    nnzGtG = (int)receivingBuffer.size();

    for (auto& ii : (*m_p_defectPerSubdomains))
      m_GtG_dim += ii;
  }

  int tmpInt[] = {nnzGtG, m_GtG_dim};

  m_p_domain->hmpi.BcastInt(&tmpInt,2,m_root);

  nnzGtG = tmpInt[0];
  m_GtG_dim = tmpInt[1];


#if DBG > 2
  std::cout << "m_GtG_dim:    " << m_GtG_dim<< '\n'; 
  std::cout << "nnzGtG:       " << nnzGtG << '\n'; 
#endif




#if !defined(GTG_ALL_NODES)
  if (m_p_domain->GetRank()  == m_root)
  {
#endif
    m_cumulativeDefectPerSubdomains.resize((*m_p_defectPerSubdomains).size(),0);
    for (int ii = 0; ii < (int) (*m_p_defectPerSubdomains).size() - 1 ; ii++)
    {
      m_cumulativeDefectPerSubdomains[ii+1] = 
        m_cumulativeDefectPerSubdomains[ii] + (*m_p_defectPerSubdomains)[ii];
    }
#if !defined(GTG_ALL_NODES)
  }
#endif

  if (m_p_domain->GetRank()  == m_root)
  {

    I_COO_GtG.resize(nnzGtG);
    J_COO_GtG.resize(nnzGtG);
    V_COO_GtG.resize(nnzGtG);



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

      offsetRow = m_cumulativeDefectPerSubdomains[diagRank];
      offsetCol = m_cumulativeDefectPerSubdomains[diagRank];

      _placeBlockInGlobalGtG(I_COO_GtG,J_COO_GtG,V_COO_GtG,
       //   trGtG, 
       GtG_local_diag,cntGtG,offsetRow ,offsetCol);



#if DBG > 4
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


        offsetRow = m_cumulativeDefectPerSubdomains[diagRank];
        offsetCol = m_cumulativeDefectPerSubdomains[offdiagRank];

        _placeBlockInGlobalGtG(I_COO_GtG,J_COO_GtG,V_COO_GtG,
           // trGtG, 
            GtG_local_offdiag,cntGtG,offsetRow, offsetCol);

#if DBG > 4
        std::cout <<  '[' << in <<',' << offdiagRank << ']' << '\n';
        std::cout << GtG_local_offdiag << '\n';
        std::cout << '\n';
#endif

        cntAllElements += defectDiag * defectOffdiag;
      }
    }
  }


#if defined(GTG_ALL_NODES)
  I_COO_GtG.resize(nnzGtG);
  J_COO_GtG.resize(nnzGtG);
  V_COO_GtG.resize(nnzGtG);

  m_p_domain->hmpi.BcastInt(I_COO_GtG.data(),nnzGtG,m_root);
  m_p_domain->hmpi.BcastInt(J_COO_GtG.data(),nnzGtG,m_root);
  m_p_domain->hmpi.BcastDbl(V_COO_GtG.data(),nnzGtG,m_root);
#endif

#if !defined(GTG_ALL_NODES)
  if (m_p_domain->GetRank()  == m_root)
  {
#endif
    std::vector<T> trGtG(nnzGtG,T(0,0,0));

    for (int iG = 0; iG < nnzGtG; iG++)
      trGtG[iG] = T(I_COO_GtG[iG],J_COO_GtG[iG],V_COO_GtG[iG]);

    m_spmatGtG.resize(m_GtG_dim, m_GtG_dim);
    m_spmatGtG.setFromTriplets(trGtG.begin(),trGtG.end());
  //
#if DBG > 3
    std::string fname = "GtG"; fname += std::to_string(m_p_domain->GetRank()) + ".txt";
    tools::printMatrix(m_spmatGtG,fname);
#endif


#if !defined(GTG_ALL_NODES)
  }
#endif


}


void InterfaceOperatorB::mult_invGtG(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  int in_cols = in.cols();
  int in_rows = in.rows();

  if ((*m_p_defectPerSubdomains)[m_p_domain->GetRank()] != in_rows)
    std::runtime_error("input vector has inccorrect number of rows");

  std::vector<int> sizePerRank(0);
  std::vector<int> sizePerRankColumPtr(0);
  std::vector<double> receivingBuffer(0);

  int nSubdomains = m_p_domain->GetNumberOfSubdomains();

  {

    sizePerRank.resize(nSubdomains,-1);
    sizePerRankColumPtr.resize(nSubdomains,-1);

    int cntAllElements(0);
    for (int i_rank = 0; i_rank < nSubdomains; i_rank++ )
    {
      int defectOnCurrentRank = (*m_p_defectPerSubdomains)[i_rank];
      sizePerRank[i_rank] = in_cols * defectOnCurrentRank;
      sizePerRankColumPtr[i_rank] = cntAllElements;
      cntAllElements += in_cols * defectOnCurrentRank;
    }

    sizePerRankColumPtr.push_back(cntAllElements);
    receivingBuffer.resize(cntAllElements,-1);
  }

  m_p_domain->hmpi.GathervDbl(in.data(), in.size(),
    receivingBuffer.data(), sizePerRank.data(),
    sizePerRankColumPtr.data(),  m_root);

  Eigen::MatrixXd sol_glb(0,0);
  if (m_p_domain->GetRank()  == m_root)
  {


    Eigen::MatrixXd rhs_glb(m_GtG_dim,in_cols);

    int offset(0);
    int cntAllElements(0);
    for (int i_rank = 0; i_rank < nSubdomains; i_rank++ )
    {
      int defectOnCurrentRank = (*m_p_defectPerSubdomains)[i_rank];

      if (defectOnCurrentRank > 0)
      {
        Eigen::Map<Eigen::MatrixXd> tmp(receivingBuffer.data() + cntAllElements,
            defectOnCurrentRank, in_cols);


        rhs_glb.block(offset, 0,
          defectOnCurrentRank, in_cols) = tmp;
      }

      offset += defectOnCurrentRank;
      cntAllElements +=  defectOnCurrentRank * in_cols;
    }

    _solve(rhs_glb,sol_glb);

  }

  out.resize(in.rows(),in.cols());
  m_p_domain->hmpi.ScattervDbl(sol_glb.data(), sizePerRank.data(),
      sizePerRankColumPtr.data(), out.data(), out.size(), m_root);
}



Eigen::MatrixXd InterfaceOperatorB::Projection(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  // Although G = (-1) * B * R
  // in step 1) and 4) ther's no multiplying by -1
  //


  Eigen::MatrixXd Bt_in, Gt_in, invGtG_Gt_in, R_invGtG_Gt_in, G_invGtG_Gt_in;
  auto &R = *(m_p_domain->GetStiffnessMatrix()->GetKernel());

  // 0) y0 = Bt * in
  multBt(in,Bt_in);

  // 1) y1 = (-1) * Rt * Bt * in = Gt * in
    Gt_in =  R.transpose() * Bt_in;  Bt_in.resize(0,0);

  // 2) y2 = alpha =  inv(GtG) * y1 = inv(GtG) * Gt * in
  mult_invGtG(Gt_in,invGtG_Gt_in);  Gt_in.resize(0,0);

  // 3) y3 = R * y2 = R * inv(GtG) * Gt * in 
  R_invGtG_Gt_in =  R * invGtG_Gt_in;    

  // 4) y4 = (-1) * B * R * y2 = G * inv(GtG) * Gt * in
  multB(R_invGtG_Gt_in,G_invGtG_Gt_in); R_invGtG_Gt_in.resize(0,0);

  // out = (I - G * inv(GtG) * Gt) * in  
  out = in - G_invGtG_Gt_in;

  return invGtG_Gt_in;

}



void InterfaceOperatorB::mult_F(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  Eigen::MatrixXd Bt_in, Kplus_Bt_in;
  auto *Kplus = m_p_domain->GetStiffnessMatrix();
  multBt(in, Bt_in);
  Kplus->solve(Bt_in,Kplus_Bt_in);   Bt_in.resize(0,0);
  multB(Kplus_Bt_in, out);  Kplus_Bt_in.resize(0,0);

}

void InterfaceOperatorB::_SetScaling()
{


#if DBG > 1
  std::cout << " --- scaling --- \n";
#endif

  auto interfaces = m_p_domain->GetInterfaces();

  int offset = 0;
  m_scaling.resize(m_p_domain->GetNumberOfDualDOFs(),0);

  for (auto& iItf : interfaces)
  {
    int nDOFs = (int) iItf.m_interfaceDOFs.size();
    for (int row = 0; row < nDOFs;row++)
    {
      m_scaling[row + offset] = (*m_p_domain->GetMultiplicity())[iItf.m_interfaceDOFs[row]];
    }
    offset += nDOFs;
  }


#if DBG > 1
  std::cout << " --- start --- \n";
  for (auto& ii : m_scaling)
    std::cout << ii << ' ';
  std::cout << '\n';
#endif

}



void InterfaceOperatorB::Scaling(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{
  out.resize(in.rows(),in.cols());

  for (int row = 0; row < in.rows(); row++)
    out.row(row) = in.row(row) / m_scaling[row];

}

void InterfaceOperatorB::Scaling(Eigen::MatrixXd& inout)
{

  for (int row = 0; row < inout.rows(); row++)
    inout.row(row) /=  m_scaling[row];

}




void InterfaceOperatorB::SFETI_Beta(Eigen::MatrixXd& inZs)
{


  // Zs projected residual from my rank

  auto kerK = m_p_domain->GetStiffnessMatrix()->GetKernel();
  auto interfaces = m_p_domain->GetInterfaces();
//
  int myDefect = m_p_domain->GetStiffnessMatrix()->GetDefect(); 
  std::cout <<  "XXX = " <<m_cumulativeDefectPerSubdomains.size() << '\n';
  int offset_myRank = m_cumulativeDefectPerSubdomains[m_p_domain->GetRank() ];
//
  // buffers
  std::vector<Eigen::MatrixXd> G_myRank(interfaces.size(),Eigen::MatrixXd(0,0));
  std::vector<Eigen::MatrixXd> G_neighb(interfaces.size(),Eigen::MatrixXd(0,0));

  int colsZs= inZs.cols();

  Eigen::MatrixXd RHS_Beta_s = Eigen::MatrixXd::Zero(m_GtG_dim,colsZs);


  int offsetDual(0);


// each rank computes whole block row 
for (int cntI = 0; cntI < (int)interfaces.size(); cntI++)
{

    Interface &iItf = interfaces[cntI];
    int nDOFsIntf = (int) iItf.m_interfaceDOFs.size();

    ////////////////////////////////////////////////////
    // CONVENTION: subdomain with higher rank has
    // operator B (or G) with negative coefficients (-1)
    ////////////////////////////////////////////////////
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int neqIntf = (int) iItf.m_interfaceDOFs.size();

    int neighbDefect = (*m_p_defectPerSubdomains)[iItf.GetNeighbRank()];
    int offset_neighb = m_cumulativeDefectPerSubdomains[iItf.GetNeighbRank()];


    G_myRank[cntI].resize(neqIntf,myDefect);
    G_neighb[cntI].resize(neqIntf,neighbDefect);

    for (int row = 0; row < neqIntf; row++)
      G_myRank[cntI].row(row) = (-1) * scale * (*kerK).row(iItf.m_interfaceDOFs[row]);

    m_p_domain->hmpi.IsendDbl(G_myRank[cntI].data(),neqIntf * myDefect, iItf.GetNeighbRank());
    m_p_domain->hmpi.IrecvDbl(G_neighb[cntI].data(),neqIntf * neighbDefect, iItf.GetNeighbRank());
    m_p_domain->hmpi.Wait();

    // myRank contribution
    RHS_Beta_s.block(offset_myRank,0,myDefect,colsZs) +=
      G_myRank[cntI].transpose() * inZs.block(offsetDual,0,nDOFsIntf,colsZs);
    // neighbr contribution
    RHS_Beta_s.block(offset_neighb,0,neighbDefect,colsZs) =
      G_neighb[cntI].transpose() * inZs.block(offsetDual,0,nDOFsIntf,colsZs);


    offsetDual += nDOFsIntf;
}

  std::cout << "RHS_Beta_s\n";
  std::cout << RHS_Beta_s;
  std::cout << "\nend\n";

  Eigen::MatrixXd Beta_s = m_pardisoSolver.solve(RHS_Beta_s);
  
  std::cout << "Beta_s\n";
  std::cout << Beta_s;
  std::cout << "\nend\n";

// empty buffers
  for (int cntI = 0; cntI < (int)interfaces.size(); cntI++){
    G_myRank[cntI].resize(0,0); G_neighb[cntI].resize(0,0);
  }


}


void InterfaceOperatorB::printInterfaceDOFs(std::string name)
{

  if (name.size()==0) name = "matB";
  std::string fname = name + std::to_string(m_p_domain->GetRank()) + ".txt";

  int nLambda = m_p_domain->GetNumberOfDualDOFs();

  Eigen::MatrixXd _matB(nLambda,4);

  auto interfaces = m_p_domain->GetInterfaces();

  int offset = 0;

  int cnt(0);

  for (auto& iItf : interfaces)
  {
    ////////////////////////////////////////////////////
    // CONVENTION: subdomain with higher rank has
    // operator B (or G) with negative coefficients (-1)
    ////////////////////////////////////////////////////
    double scale(1.0);
    if (m_p_domain->GetRank() > iItf.GetNeighbRank())
      scale *= -1;

    int nDOFs = (int) iItf.m_interfaceDOFs.size();
    for (int row = 0; row < nDOFs;row++)
    {
      _matB(cnt,0) = iItf.GetNeighbRank();
      _matB(cnt,1) = iItf.m_interfaceDOFs[row];
      _matB(cnt,2) = scale;
      _matB(cnt,3) = m_scaling[cnt];
      cnt++;
    }
    offset += nDOFs;
  }

  tools::printMatrix(_matB,fname);

}





void InterfaceOperatorB::printNeighboursRanks(std::string name)
{

  if (name.size()==0) name = "neighRank";
  std::string fname = name + std::to_string(m_p_domain->GetRank()) + ".txt";

  auto interfaces = m_p_domain->GetInterfaces();

  std::vector<int> nghb(0);
  for (auto& iItf : interfaces)
  {
    nghb.push_back(iItf.GetNeighbRank());
  }

  tools::printArray(nghb,fname);

}




void InterfaceOperatorB::m_releaseBuffers()
{
  for(auto& ii : m_sendBuffers) ii.resize(0,0);
}

