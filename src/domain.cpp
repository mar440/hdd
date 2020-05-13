#include "../include/domain.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/stiffnessMatrix.hpp"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <set>
#include <map>

#include <chrono>


#define TAG0 100
#define TAG1 101

Domain::Domain(MPI_Comm* _pcomm): hmpi(_pcomm)
{
  m_comm =  *_pcomm;
  m_Init();
}

Domain::~Domain()
{

  if (m_p_stiffnessMatrix) delete m_p_stiffnessMatrix;

}

void Domain::m_Init()
{
  m_multiplicity.resize(0);
  m_neighboursRanks.resize(0);
  m_l2g.resize(0);
  m_neqPrimal = 0;
  m_neqDual= 0;
  m_DirichletDOFs.resize(0);

  m_I_DirichletPrecondDOFs.resize(0);
  m_B_DirichletPrecondDOFs.resize(0);

  MPI_Comm_rank(m_comm,&m_mpirank);
  MPI_Comm_size(m_comm,&m_mpisize);

}

void Domain::InitStiffnessMatrix(PRECONDITIONER_TYPE precondType)
{
  m_p_stiffnessMatrix = new StiffnessMatrix(&m_g2l,m_mpirank, precondType);
}

void Domain::SetMappingLoc2Glb(std::vector<int>& _l2g)
{
  m_l2g = _l2g;
  m_neqPrimal = (int) m_l2g.size();
  m_SetMappingGlb2Loc();
}

void Domain::m_SetMappingGlb2Loc()
{
  // maping DOF: global to local
  for (int i = 0; i < m_neqPrimal; i++)
    m_g2l[m_l2g[i]] = i;
}

void Domain::SetNeighboursRanks(const std::vector<int>& ranks) 
{

  m_neighboursRanks = ranks;
  
#if DBG > 0
  m_dbg_printNeighboursRanks();
#endif

}


void Domain::SetDirichletDOFs(std::vector<int>& glbDirDOFs)
{

  std::cout << "Dirichlet: glb -> loc \n";
  std::cout <<"m_g2l.size(): " << m_g2l.size() << '\n';

  for (auto& id : glbDirDOFs) 
  {
    auto it = m_g2l.find(id);
    if (  it != m_g2l.end())
      m_DirichletDOFs.push_back(it->second);
  }

  std::cout << "number of Dir. DOFs: " << m_DirichletDOFs.size() << '\n';
}




void Domain::m_dbg_printNeighboursRanks()
{

  std::cout << "neighboursRanks(Domain): ";
  for (auto& indx : m_neighboursRanks)
    std::cout << indx << ' ';
  std::cout << '\n';
}


void Domain::_SearchNeighbours(std::vector<int>& inout)
{


  // TODO make it parallel - independently on domain decomposition
  HDDTRACES

  auto startTime = std::chrono::steady_clock::now();

  HDDTRACES
  // get max DOF index in global numbering
  auto maxDofIndOnDomain =
    max_element(m_l2g.begin(), m_l2g.end());

  int maxDofInd =(int) *maxDofIndOnDomain + 1;
  hmpi.GlobalInt(&maxDofInd,1,MPI_MAX);
  if (m_mpirank == 0) std::cout << "maxDofInd: " << maxDofInd<< '\n';

  HDDTRACES
  std::vector<int> multiplicity(maxDofInd,0);

  HDDTRACES
  for (auto& iii : m_l2g)
    multiplicity[iii] = 1;

//  std::cout<< "MPC_SIZE: " << multiplicity.size() << '\n';
//  std::cout << "SSS=====\n";
//  for (auto& iii : multiplicity) std::cout << ' ' << iii;
//  std::cout << "\nEEE=====\n";

  HDDTRACES

//hmpi.GlobalInt(multiplicity,MPI_SUM);
  hmpi.GlobalInt(multiplicity.data(),multiplicity.size(),MPI_SUM);

  HDDTRACES


  std::vector<int> dofOnIntf(0);
  for (auto& iii : m_l2g)
  {
    if (multiplicity[iii] > 1)
      dofOnIntf.push_back(iii);
  }
//  std::cout << "SSS=====\n";
//  for (auto& iii : dofOnIntf) std::cout << iii << '\n';
//  std::cout << "EEE=====\n";

  HDDTRACES
  std::vector<int> numberOfDofOnIntfZ0(m_mpisize,0);
  std::vector<int> ptrZ0(m_mpisize,0);

  int numberOfDofOnIntf = (int) dofOnIntf.size();
  numberOfDofOnIntfZ0[m_mpirank] = numberOfDofOnIntf;

  int root(0);
  hmpi.GatherInt(&numberOfDofOnIntf, 1,
               numberOfDofOnIntfZ0.data(), 1 , root);

  HDDTRACES
  int sumIntfDofs(0);
  if (m_mpirank == root)
  {
//    std::cout << "SSS=====\n";
    for (auto& iii : numberOfDofOnIntfZ0)
    {
  //    std::cout << iii << '\n';
      sumIntfDofs +=  iii;
    }
//    std::cout << "EEE=====\n";
  }

  HDDTRACES
// ----------------------------------
  int edgeSum=0;

  std::vector<int>  intfDofsOnRoot;
  if (m_mpirank  == root){
    intfDofsOnRoot.resize(sumIntfDofs,-1);

    ptrZ0.resize(m_mpisize + 1);
    for(int i=0; i < m_mpisize; ++i) {
        ptrZ0[i]=edgeSum;
        edgeSum+=numberOfDofOnIntfZ0[i];
    }
    ptrZ0[m_mpisize] = edgeSum;
  }
// ----------------------------------  
  HDDTRACES

  hmpi.GathervInt(dofOnIntf.data(), dofOnIntf.size(), // to send
      intfDofsOnRoot.data(),
      numberOfDofOnIntfZ0.data(),ptrZ0.data(),
      root);

  std::vector<std::set<int>> neighbours;
  std::map<int,std::list<int>> dofAndDomain;
  int sumOfInterfacesGlobal(0);

  HDDTRACES
  if (m_mpirank  == root)
  {
    for (int isub = 0; isub < m_mpisize; isub++)
    {
      for (int jj = ptrZ0[isub]; jj < ptrZ0[isub+1] ; jj++)
      {
//        std::cout <<  intfDofsOnRoot[jj] << ' ';
        dofAndDomain[intfDofsOnRoot[jj]].push_back(isub);
      }
 //     std::cout << '\n';
    }

    //for (auto& imap : dofAndDomain)
    //{
    //  std::cout << imap.first << ": ";
    //  for (auto& ilist : imap.second)
    //    std::cout << ilist << ' ';
    //  std::cout << '\n';
    //}

    m_listOfNeighboursColumPtr.resize(m_mpisize+1,0);
    neighbours.resize(m_mpisize);
    for (int isub = 0; isub < m_mpisize; isub++)
    {
      for (int jj = ptrZ0[isub]; jj < ptrZ0[isub+1] ; jj++)
      {
        for (int& kk : dofAndDomain[intfDofsOnRoot[jj]])
        {
          neighbours[isub].insert(kk);
//          std::cout << isub << ":" << intfDofsOnRoot[jj]<<":" << kk << "\n";
        }
      }
      auto it = neighbours[isub].find(isub);
      neighbours[isub].erase(it);
      sumOfInterfacesGlobal += neighbours[isub].size();
      m_listOfNeighboursColumPtr[isub + 1] = sumOfInterfacesGlobal;
    }

    std::cout << "neighb....\n";
    for (int isub = 0; isub < m_mpisize; isub++)
    {
      std::cout << isub << ": ";
      for (auto& ns : neighbours[isub])
      {
        std::cout << ns << ' ';
      }
      std::cout << '\n';
    }
  } // rank == root 

  HDDTRACES
  std::vector<int> numberOfNeighboursRoot;
  int numberOfNeighboursLocal(0);

  m_listOfNeighbours.resize(sumOfInterfacesGlobal);

  if (m_mpirank  == root)
  {
    int cnt(0);
    numberOfNeighboursRoot.resize(m_mpisize,-1);
    for (int isub = 0; isub < m_mpisize; isub++)
    {
      numberOfNeighboursRoot[isub] = neighbours[isub].size();
      for (auto &kface : neighbours[isub])
        m_listOfNeighbours[cnt++]  = kface;
    }
  }

  hmpi.ScatterInt(numberOfNeighboursRoot.data(),1, 
      &numberOfNeighboursLocal, 1, root);

  std::cout << "number of neighbours is: " << numberOfNeighboursLocal << '\n';
  std::vector<int> listOfInterfacesLocal(numberOfNeighboursLocal,-1);

  std::cout << "------------\n";
  for (auto& intf : m_listOfNeighboursColumPtr)
    std::cout << intf << ' ';
  std::cout << '\n';

  HDDTRACES
//
  hmpi.ScattervInt(m_listOfNeighbours.data(), numberOfNeighboursRoot.data(),
      m_listOfNeighboursColumPtr.data(), 
      listOfInterfacesLocal.data(),listOfInterfacesLocal.size(),root);

  for (auto& intf : listOfInterfacesLocal)
    std::cout << intf << ' ';
  std::cout << '\n';

  inout = listOfInterfacesLocal;

  auto endTime = std::chrono::steady_clock::now();
  std::chrono::duration<double> elapsed_seconds = endTime-startTime;
  std::cout << std::fixed << std::setprecision(2) << 
      "Search neighbours: " << elapsed_seconds.count() << " s\n";
  HDDTRACES
}

void Domain::SetInterfaces()
{

  HDDTRACES
  // if "neighbours" ranks not known apriori
  if ((int)m_neighboursRanks.size() == 0)
    _SearchNeighbours(m_neighboursRanks);

  HDDTRACES

  int nInterf = (int)m_neighboursRanks.size();

  m_interfaces.resize(nInterf);
  for (int iR = 0; iR < nInterf; iR++)
  {
    m_interfaces[iR].SetNeighbRank(m_neighboursRanks[iR]);
    m_intf_g2l[m_neighboursRanks[iR]] = iR;
  }
  HDDTRACES
  {

    MPI_Request requests[nInterf];
    MPI_Status statuses[nInterf];
    HDDTRACES
    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {

      Interface& i_intfc = m_interfaces[intfId];
      HDDTRACES

      int neighRank = i_intfc.GetNeighbRank();

      MPI_Isend(m_l2g.data(),m_l2g.size(),MPI_INT,neighRank,TAG0,
              hmpi.GetComm(),&requests[intfId]);
      HDDTRACES
    }

    HDDTRACES
    int msg_avail = -1;
    MPI_Status status1, status2;
    int neq_neighb(0);
    int cntRecvMsg(0);

    while(true)
    {

      MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, hmpi.GetComm(),&msg_avail, &status1);

      HDDTRACES
      if (msg_avail!=0)
      {
        HDDTRACES
        MPI_Get_count(&status1, MPI_INT,&neq_neighb);
        int neighRank = status1.MPI_SOURCE;

        std::vector<int>  neighb_l2g(neq_neighb,-1);
        MPI_Recv(neighb_l2g.data(),neq_neighb,MPI_INT, neighRank,TAG0,
            hmpi.GetComm(), &status2);

        std::vector<int> intersection = 
          tools::intersection(m_l2g, neighb_l2g);

        int neqIntf = intersection.size();

        int intfId = m_intf_g2l[neighRank];

        Interface& i_intfc = m_interfaces[intfId];
        i_intfc.SetNeighbNumbOfEqv(neq_neighb);

        i_intfc.m_interfaceDOFs.resize(neqIntf);
        for (int dof = 0; dof < neqIntf; dof ++)
          i_intfc.m_interfaceDOFs[dof] =  m_g2l[intersection[dof]];

        std::cout << "neighRank:" << neighRank << ", neq: " <<  neq_neighb <<std::endl;

        HDDTRACES
        if (++cntRecvMsg == nInterf) break;
        HDDTRACES
      }
    }
    MPI_Wait(requests,statuses);
  }

  HDDTRACES

  // set weight and get dual number of dofs
  m_neqDual = 0;
  m_multiplicity.resize(m_neqPrimal,1);

  for (auto& i_intfc : m_interfaces)
  {
    for(auto& ddof : i_intfc.m_interfaceDOFs)
      m_multiplicity[ddof]++;

    m_neqDual += i_intfc.m_interfaceDOFs.size();

  }

  for (auto& im :  m_multiplicity) im = sqrt(im);


  // set offset
//  for (auto& intf : m_interfaces)
  

#if DBG > 3
  for (auto& im :  m_multiplicity)
    std::cout << im << ' ';
  std::cout << std::endl;

#endif

}

void Domain::_SetInterfaces()
{

  // LEGACY ...

  // SET m_interfaces done in 3 steps
  //    1 - number of DOFs exchanging between neigbhours
  //    2 - sending DOFs from lower to higher rank and
  //        determine the intersection
  //    3 - sending intersection (from higher to lower rank)

  HDDTRACES

  // if "neighbours" ranks not known apriori
  if ((int)m_neighboursRanks.size() == 0)
    _SearchNeighbours(m_neighboursRanks);

  HDDTRACES

  int nInterf = (int)m_neighboursRanks.size();

  m_interfaces.resize(nInterf);
  for (int iR = 0; iR < nInterf; iR++)
  {
    m_interfaces[iR].SetNeighbRank(m_neighboursRanks[iR]);
    m_intf_g2l[m_neighboursRanks[iR]] = iR;
  }

    

  HDDTRACES
  
  //##################################
  // 1 get my neighbours neq (number of DOFs) 
  {

    std::vector<int> neq_neighbs(nInterf,0);
    std::vector<MPI_Request> requests(0);

    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {
      Interface& i_intfc = m_interfaces[intfId];

      int neighRank = i_intfc.GetNeighbRank();

      MPI_Send(&m_neqPrimal,1,MPI_INT,neighRank,TAG0,
          hmpi.GetComm());


    }
//    for (int intfId = 0 ; intfId < nInterf; intfId++ )

    int flag = -1;
    MPI_Request request;
    MPI_Status status;
    int neq_neighb(0);

    int cntRecvMsg(0);

    while(true)
    {
      if (flag!=0)
      {
        MPI_Irecv(&neq_neighb,1,MPI_INT,MPI_ANY_SOURCE,MPI_ANY_TAG,
            hmpi.GetComm(), &request);
        flag = 0;
      }
      MPI_Test(&request,&flag,&status);
      if (flag!=0)
      {
        int neighRank = status.MPI_SOURCE;
        int intfId = m_intf_g2l[neighRank];
        neq_neighbs[intfId] = neq_neighb;
        Interface& i_intfc = m_interfaces[intfId];
        i_intfc.SetNeighbNumbOfEqv(neq_neighbs[intfId]);

        std::cout << "neighRank:" << neighRank << ", neq: " <<  neq_neighb <<std::endl;

        if (++cntRecvMsg == nInterf) break;
      }

    }


    hmpi.Barrier();

//    for (int intfId = 0 ; intfId < nInterf; intfId++ )
//    {
//      Interface& i_intfc = m_interfaces[intfId];
//      i_intfc.SetNeighbNumbOfEqv(neq_neighbs[intfId]);
//    }
  }


  //##################################
  // 2 get neigbour's global DOFs to
  // determine common DOFs (intersection)
  std::vector<std::vector<int>> 
    intersections(nInterf,std::vector<int>(0));

  {

    std::vector<MPI_Request> requests(0);

    std::vector<std::vector<int>> 
      neighbs_l2g(nInterf,std::vector<int>(0));

    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {

      Interface& i_intfc = m_interfaces[intfId];

      int neq_neighb  = i_intfc.GetNeighbNumbOfEqv();
      int neighRank   = i_intfc.GetNeighbRank();

      MPI_Request newRequest;

      if (m_mpirank > neighRank)
      {
        neighbs_l2g[intfId].resize(neq_neighb);
        MPI_Irecv(neighbs_l2g[intfId].data(),neq_neighb,MPI_INT, neighRank,TAG0,
            hmpi.GetComm(),&newRequest);
      }
      else
      {
        MPI_Isend(m_l2g.data(),m_l2g.size(),MPI_INT,neighRank,TAG0,
            hmpi.GetComm(), &newRequest);
      }
      requests.push_back(newRequest);
    }

    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall((int)requests.size(),requests.data(),status.data());

    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {
      Interface& i_intfc = m_interfaces[intfId];
      int neighRank   = i_intfc.GetNeighbRank();

      if (m_mpirank > neighRank)
      {
        intersections[intfId] = tools::intersection(m_l2g, neighbs_l2g[intfId]);
        neighbs_l2g[intfId].resize(0);
        neighbs_l2g[intfId].shrink_to_fit();
      }
    }

  }

//  hmpi.Barrier();

  //##################################
  // 3 send intersection DOFs to lower rank 
  std::vector<int> neqInterfaces(nInterf,0);
  {

    std::vector<MPI_Request> requests(0);
    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {

      Interface& i_intfc = m_interfaces[intfId];
      int neighRank   = i_intfc.GetNeighbRank();

      MPI_Request newRequest;
      if (m_mpirank > neighRank)
      {
        neqInterfaces[intfId] = intersections[intfId].size();
        MPI_Isend(&(neqInterfaces[intfId]),1,MPI_INT,neighRank,TAG0,
            hmpi.GetComm(), &newRequest);
      }
      else
      {
        MPI_Irecv(&(neqInterfaces[intfId]),1,MPI_INT, neighRank,TAG0,
            hmpi.GetComm(),&newRequest);
      }
      requests.push_back(newRequest);
    }

    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall((int)requests.size(),requests.data(),status.data());
  }

  //##################################
  // 3 send intersection DOFs to lower rank 
  {

    std::vector<MPI_Request> requests(0);

    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {

      Interface& i_intfc = m_interfaces[intfId];
      int neighRank   = i_intfc.GetNeighbRank();

      std::cout << "neighbRank: " << i_intfc.GetNeighbRank();
      std::cout << " -- neqInterface = " << neqInterfaces[intfId] << std::endl;

      MPI_Request newRequest;

      if (m_mpirank > neighRank)
      {
        MPI_Isend(intersections[intfId].data(),neqInterfaces[intfId],MPI_INT,
            neighRank,TAG0, hmpi.GetComm(), &newRequest);
      }
      else
      {
        intersections[intfId].resize(neqInterfaces[intfId],-1);
        MPI_Irecv(intersections[intfId].data(),neqInterfaces[intfId],MPI_INT, 
            neighRank,TAG0, hmpi.GetComm(),&newRequest);
      }

      requests.push_back(newRequest);

    }


    std::vector<MPI_Status> status(requests.size());
    MPI_Waitall((int)requests.size(),requests.data(),status.data());

    for (int intfId = 0 ; intfId < nInterf; intfId++ )
    {
      Interface& i_intfc = m_interfaces[intfId];
      i_intfc.m_interfaceDOFs.resize(neqInterfaces[intfId]);
      for (int dof = 0; dof < neqInterfaces[intfId]; dof ++)
        i_intfc.m_interfaceDOFs[dof] =  m_g2l[intersections[intfId][dof]];

      intersections[intfId].resize(0);
      intersections[intfId].shrink_to_fit();

    }

  }

  HDDTRACES

  // set weight and get dual number of dofs
  m_neqDual = 0;
  m_multiplicity.resize(m_neqPrimal,1);

  for (auto& i_intfc : m_interfaces)
  {
    for(auto& ddof : i_intfc.m_interfaceDOFs)
      m_multiplicity[ddof]++;

    m_neqDual += i_intfc.m_interfaceDOFs.size();

  }

  for (auto& im :  m_multiplicity) im = sqrt(im);

#if DBG > 3
  for (auto& im :  m_multiplicity)
    std::cout << im << ' ';
  std::cout << std::endl;

#endif

}
void Domain::_SetDirichletPrecondDOFs()
{
  std::vector<bool> isDualDOF(m_neqPrimal,false);
  for (auto& i_intfc : m_interfaces)
  {
    for(auto& ddof : i_intfc.m_interfaceDOFs)
    {
      isDualDOF[ddof] = true;
    }
  }


  //for (auto idof : isDualDOF)
  for (int id = 0; id < (int)isDualDOF.size(); id++)
  {
    if (isDualDOF[id])
      m_B_DirichletPrecondDOFs.push_back(id);
    else
      m_I_DirichletPrecondDOFs.push_back(id);

  }
}

void Domain::HandlePreconditioning()
{
    
  _SetDirichletPrecondDOFs();
 
  m_p_stiffnessMatrix->SetDirichletPrecond(
        m_I_DirichletPrecondDOFs, 
        m_B_DirichletPrecondDOFs);

}


/////////////////////////////////////////////////////////////////////////
// DBG  DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG DBG
/////////////////////////////////////////////////////////////////////////

void Domain::m_dbg_print_l2g()
{
  std::cout << "l2g(" << m_l2g.size() << ")"<< '\n';
  int cnt(0);
  for (auto& id : m_l2g)
  {
    std::cout << id << ' ';
    if ((++cnt)%10 == 0)
      std::cout <<'\n';
  }
  std::cout <<'\n';
}



