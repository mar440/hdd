#include "../include/domain.hpp"
#include "../include/linearAlgebra.hpp"
#include "../include/stiffnessMatrix.hpp"

#include <iostream>
#include <iomanip>
#include <math.h>
#include <set>
#include <map>

#include <chrono>

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

  int cntD(0);

  for (auto& id : glbDirDOFs) 
  {

    auto it = m_g2l.find(id);

    if (  it != m_g2l.end())
    {
      m_DirichletDOFs.push_back(it->second); 
      cntD++;
    }
  }

  std::cout << "number of Dir. DOFs: " << cntD << '\n';
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

  auto startTime = std::chrono::steady_clock::now();

  // get max DOF index in global numbering
  auto maxDofIndOnDomain =
    max_element(m_l2g.begin(), m_l2g.end());

  int maxDofInd = *maxDofIndOnDomain;
  hmpi.GlobalInt(&maxDofInd,1,MPI_MAX);
  if (m_mpirank == 0) std::cout << "maxDofInd: " << maxDofInd<< '\n';

  std::vector<int> multiplicity(maxDofInd,0);

  for (auto& iii : m_l2g)
    multiplicity[iii] = 1;

  hmpi.GlobalInt(multiplicity.data(),multiplicity.size(),MPI_SUM);

//  std::cout << "SSS=====\n";
//  for (auto& iii : multiplicity)
//    std::cout << iii << '\n';
//  std::cout << "EEE=====\n";

  std::vector<int> dofOnIntf(0);
  for (auto& iii : m_l2g)
  {
    if (multiplicity[iii] > 1)
      dofOnIntf.push_back(iii);
  }
//  std::cout << "SSS=====\n";
//  for (auto& iii : dofOnIntf) std::cout << iii << '\n';
//  std::cout << "EEE=====\n";

  std::vector<int> numberOfDofOnIntfZ0(m_mpisize,0);
  std::vector<int> ptrZ0(m_mpisize,0);

  int numberOfDofOnIntf = (int) dofOnIntf.size();
  numberOfDofOnIntfZ0[m_mpirank] = numberOfDofOnIntf;

  int root(0);
  hmpi.GatherInt(&numberOfDofOnIntf, 1,
               numberOfDofOnIntfZ0.data(), 1 , root);

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

  hmpi.GathervInt(dofOnIntf.data(), dofOnIntf.size(), // to send
      intfDofsOnRoot.data(),
      numberOfDofOnIntfZ0.data(),ptrZ0.data(),
      root);

  std::vector<std::set<int>> neighbours;
  std::map<int,std::list<int>> dofAndDomain;
  int sumOfInterfacesGlobal(0);

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


}

void Domain::SetInterfaces()
{



  if ((int)m_neighboursRanks.size() == 0)
  {
    _SearchNeighbours(m_neighboursRanks);
  }

  int nInterf = (int)m_neighboursRanks.size();

  m_interfaces.resize(nInterf);
  for (int iR = 0; iR < nInterf; iR++)
    m_interfaces[iR].SetNeighbRank(m_neighboursRanks[iR]);



  for (int i_ = 0 ; i_ < nInterf; i_++ )
  {

    Interface& i_intfc = m_interfaces[i_];

    int neighRank = i_intfc.GetNeighbRank();
    int neq_neighb(0);


    hmpi.SendInt(&m_neqPrimal,1 , neighRank);
    hmpi.RecvInt(&neq_neighb,1 , neighRank);

    i_intfc.SetNeighbNumbOfEqv(neq_neighb);
    std::vector<int> intersection(0);

    if (m_mpirank > neighRank)
    {
      // greather rank manages intersection
      std::vector<int> neighb_l2g(neq_neighb,0);
      hmpi.RecvInt(neighb_l2g.data(),neq_neighb , neighRank);
      intersection = tools::intersection(m_l2g, neighb_l2g);
    }
    else
    {
      // send my l2g to neighbour with < rank
      hmpi.SendInt(m_l2g.data(),m_l2g.size(), neighRank);
    }

    int neqInterface(0);

    if (m_mpirank > neighRank)
    {
      neqInterface = intersection.size();
      hmpi.SendInt(&neqInterface,1, neighRank);
    }
    else
    {
      hmpi.RecvInt(&neqInterface,1, neighRank);
    }

    std::cout << "neighbRank: " << neighRank;
    std::cout << " -- neqInterface = " << neqInterface << std::endl;


    if (m_mpirank > neighRank)
    {
      hmpi.SendInt(intersection.data(),neqInterface, neighRank);
    }
    else
    {
      intersection.resize(neqInterface,-1);
      hmpi.RecvInt(intersection.data(),neqInterface, neighRank);
    }

    i_intfc.m_interfaceDOFs.resize(neqInterface);
    for (int dof = 0; dof < neqInterface; dof ++)
      i_intfc.m_interfaceDOFs[dof] =  m_g2l[intersection[dof]];
//
//
//
#if DBG > 2
    std::cout  << "#interf. dofs:  ";
    for (auto& ii : i_intfc.m_interfaceDOFs)
      std::cout << ii << ' ' ;
    std::cout<< std::endl;
#endif

  }


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



