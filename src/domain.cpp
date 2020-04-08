#include "../include/domain.hpp"
#include <iostream>
#include <math.h>

#include "../include/linearAlgebra.hpp"
#include "../include/stiffnessMatrix.hpp"

Domain::Domain(MPI_Comm* _pcomm): hmpi(_pcomm)
{

  m_pcomm =  _pcomm;
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
  MPI_Comm_rank(*m_pcomm,&m_mpirank);
  MPI_Comm_size(*m_pcomm,&m_mpisize);

}

void Domain::InitStiffnessMatrix()
{

  m_p_stiffnessMatrix = new StiffnessMatrix(&m_g2l,m_mpirank);

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


void Domain::SetInterfaces()
{

  // m_multiplicity

  int nInterf = m_neighboursRanks.size();
  if (nInterf == 0)
    std::runtime_error("undecomposed?");
  //TODO or Dirichlet


  // allocate
  m_interfaces.resize(nInterf);
  for (int iR = 0; iR < nInterf; iR++)
    m_interfaces[iR].SetNeighbRank(m_neighboursRanks[iR]);



  for (auto& i_intfc : m_interfaces)
  {
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

    std::cout << "neqInterface= " << neqInterface << std::endl;




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

  for (auto& im :  m_multiplicity)
    im = (im - 1) * im * 0.5;
#if DBG > 2
  for (auto& im :  m_multiplicity)
    std::cout << im << ' ';
  std::cout << std::endl;

#endif



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



