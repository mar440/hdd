#include "../include/interfaceOperatorB.hpp"
#include "../include/domain.hpp"
#include <iostream>



#include "../include/stiffnessMatrix.hpp"
#include "../include/linearAlgebra.hpp"
#include <cmath>




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
  m_GtG_rows = 0;
  m_dbufToSend.resize(0,0);
  m_dbufToRecv.resize(0,0);


//  std::cout << "my rank is: " << _dom->GetRank() << std::endl;
//  std::cout << " and my interfaces ranks are: " << std::endl;
//  for (auto& ii : interfaces)
//  {
//    std::cout<<  " " << ii.GetNeighbRank() << std::endl;
//  }
}




void InterfaceOperatorB::multBt(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
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


void InterfaceOperatorB::multB(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
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

    m_dbufToRecv.resize(nDOFs , nRHS);
    m_dbufToSend.resize(nDOFs , nRHS);

    for (int row = 0; row < nDOFs; row++)
    {
      out.row(row + offset) = scale * in.row(iItf.m_interfaceDOFs[row]);
      m_dbufToSend.row(row) = scale * in.row(iItf.m_interfaceDOFs[row]);
    }

    m_p_domain->hmpi.SendDbl(m_dbufToSend.data(),nDOFs * nRHS, iItf.GetNeighbRank());
    m_p_domain->hmpi.RecvDbl(m_dbufToRecv.data(),nDOFs * nRHS, iItf.GetNeighbRank());

    for (int row = 0; row < nDOFs; row++)
      out.row(row + offset) += m_dbufToRecv.row(row);

    offset += nDOFs;
  }
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

  m_p_defectPerSubdomains =  &defectPerSubdomains;
  _FetiCoarseSpaceAssembling();

  m_pardisoSolver.analyzePattern(m_spmatGtG);
  m_pardisoSolver.factorize(m_spmatGtG);

}

void InterfaceOperatorB::solve(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  out = m_pardisoSolver.solve(in);

}




void InterfaceOperatorB::_FetiCoarseSpaceAssembling()
{

  auto kerK = m_p_domain->GetStiffnessMatrix()->GetKernel();
  int myDefect = m_p_domain->GetStiffnessMatrix()->GetDefect(); 
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

    // !!! REMAINDER in ''Gt*lambda = e'' negative sign: Gt = -Rt*Bt, e = -Rt*f
    for (int row = 0; row < neqIntf; row++)
      BR_myRank[cntI].row(row) = (-1) * scale * (*kerK).row(iItf.m_interfaceDOFs[row]);

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
  std::cout << "numberOfNeighboursRoot: << "  << numberOfNeighboursRoot.size() << '\n';

  m_p_domain->hmpi.GatherInt(&numberOfInterfacesPerRank, 1,
               numberOfNeighboursRoot.data(), 1 ,m_root);


  if (m_p_domain->GetRank() == 0)
  for (auto& ii : numberOfNeighboursRoot) std::cout << ii << ' ';
  std::cout << '\n';

  int sumOfInterfacesGlobal(0);
  m_p_domain->hmpi.ReduceInt(&numberOfInterfacesPerRank,&sumOfInterfacesGlobal,
      1, MPI_SUM, m_root);

  std::cout <<  "numberOfInterfacesPerRank: " << numberOfInterfacesPerRank<< std::endl;
  if (m_p_domain->GetRank()  == m_root)
    std::cout <<  "sumOfInterfacesGlobal: " << sumOfInterfacesGlobal << std::endl;



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



  m_p_domain->hmpi.GathervInt(
      listOfInterfacesLocal.data(), listOfInterfacesLocal.size(), // to send
      m_listOfNeighbours.data(), numberOfNeighboursRoot.data(),   // info
      m_listOfNeighboursColumPtr.data(),  m_root);


  //if (m_p_domain->GetRank()  == m_root)
  //{
  //  for (int in = 0; in < (int) numberOfNeighboursRoot.size(); in++ )
  //  {
  //    std::cout << "(" << in << "): " ;
  //    for (int jn = m_listOfNeighboursColumPtr[in];
  //        jn < m_listOfNeighboursColumPtr[in+1]; jn++)
  //    {
  //      std::cout << m_listOfNeighbours[jn]<< ' ';
  //    }
  //    std::cout << '\n';
  //  }
  //  std::cout <<  "===========================+\n";
  //}


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

  std::cout << "receivingBuffer.size() = " << receivingBuffer.size() << '\n';



  std::cout << "sendingBuffer.size(): " << sendingBuffer.size() << '\n';

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
      m_GtG_rows += ii;
  }

  int tmpInt[] = {nnzGtG, m_GtG_rows};

  m_p_domain->hmpi.BcastInt(&tmpInt,2,m_root);

  nnzGtG = tmpInt[0];
  m_GtG_rows = tmpInt[1];


  std::cout << "m_GtG_rows:   " << m_GtG_rows << '\n'; 
  std::cout << "nnzGtG:       " << nnzGtG << '\n'; 


  if (m_p_domain->GetRank()  == m_root)
  {

    I_COO_GtG.resize(nnzGtG);
    J_COO_GtG.resize(nnzGtG);
    V_COO_GtG.resize(nnzGtG);


    m_cumulativeDefectPerSubdomains.resize((*m_p_defectPerSubdomains).size(),0);
    for (int ii = 0; ii < (int) (*m_p_defectPerSubdomains).size() - 1 ; ii++)
    {
      m_cumulativeDefectPerSubdomains[ii+1] = 
        m_cumulativeDefectPerSubdomains[ii] + (*m_p_defectPerSubdomains)[ii];
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

  I_COO_GtG.resize(nnzGtG);
  J_COO_GtG.resize(nnzGtG);
  V_COO_GtG.resize(nnzGtG);
//
  m_p_domain->hmpi.BcastInt(I_COO_GtG.data(),nnzGtG,m_root);
  m_p_domain->hmpi.BcastInt(J_COO_GtG.data(),nnzGtG,m_root);
  m_p_domain->hmpi.BcastDbl(V_COO_GtG.data(),nnzGtG,m_root);

  std::vector<T> trGtG(nnzGtG,T(0,0,0));

  for (int iG = 0; iG < nnzGtG; iG++)
    trGtG[iG] = T(I_COO_GtG[iG],J_COO_GtG[iG],V_COO_GtG[iG]);

  m_spmatGtG.resize(m_GtG_rows, m_GtG_rows);
  m_spmatGtG.setFromTriplets(trGtG.begin(),trGtG.end());
//
  std::string fname = "GtG"; fname += std::to_string(m_p_domain->GetRank()) + ".txt";
  tools::printMatrix(m_spmatGtG,fname);


}

//void InterfaceOperatorB::GinvGtG(const Eigen::MatrixXd& in_loc, Eigen::MatrixXd& out_glb)
//{
//
//  // collect local Rt*X  to global Rt * X on root
//  Eigen::MatrixXd in_glb;
//  _GathervDblCoarse(in_loc,in_glb);
//  Eigen::MatrixXd invGtG_in_glb;
//
//  if (m_p_domain->GetRank() == m_root)
//  {
//    solve(in_glb,invGtG_in_glb);
//  }
//
//
//}

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


    Eigen::MatrixXd rhs_glb(m_GtG_rows,in_cols);

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

    solve(rhs_glb,sol_glb);

  }

  out.resize(in.rows(),in.cols());
  m_p_domain->hmpi.ScattervDbl(sol_glb.data(), sizePerRank.data(),
      sizePerRankColumPtr.data(), out.data(), out.size(), m_root);
}



void InterfaceOperatorB::Projection(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  Eigen::MatrixXd Bt_in, Gt_in, invGtG_Gt_in, R_invGtG_Gt_in, G_invGtG_Gt_in;

  // 0) y0 = Bt * in
  multBt(in,Bt_in);

  // 1) y1 = Rt * Bt * in = Gt * in
  auto &R = *(m_p_domain->GetStiffnessMatrix()->GetKernel());
  Gt_in =  R.transpose() * Bt_in;  Bt_in.resize(0,0); // (-1)  not multiplied  here

  // 2) y2 = inv(GtG) * y1 = inv(GtG) * Gt * in
  mult_invGtG(Gt_in,invGtG_Gt_in);  Gt_in.resize(0,0);

  // 3) y3 = R * y2 = R * inv(GtG) * Gt * in 
  R_invGtG_Gt_in =  R * invGtG_Gt_in;  invGtG_Gt_in.resize(0,0);  // (-1)  not multiplied  here

  // 4) y4 = B * R * y2 = G * inv(GtG) * Gt * in 
  multB(R_invGtG_Gt_in,G_invGtG_Gt_in); R_invGtG_Gt_in.resize(0,0);      // (-1) not multiplied here

  // out = (I - G * inv(GtG) * Gt) * in  
  out = in - G_invGtG_Gt_in;

}



void InterfaceOperatorB::mult_F(const Eigen::MatrixXd& in, Eigen::MatrixXd& out)
{

  Eigen::MatrixXd Bt_in, Kplus_Bt_in;
  
  auto *Kplus = m_p_domain->GetStiffnessMatrix();

//  multBt(in, Bt_in);
//  Kplus->solve(Bt_in,Kplus_Bt_in);   Bt_in.resize(0,0);
//  multB(Kplus_Bt_in, out);  Kplus_Bt_in.resize(0,0);

}


