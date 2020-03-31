#include "../include/domain.hpp"
#include <iostream>
#include <math.h>

#include "../include/linearAlgebra.hpp"

Domain::Domain(MPI_Comm* _pcomm): hmpi(_pcomm){

  m_pcomm =  _pcomm;
  m_Init();
}

void Domain::m_Init()
{
  m_neighboursRanks.resize(0);
  m_stiffnessMatrix.resize(0,0);
  m_l2g.resize(0);
  m_trK.resize(0);
  m_neq = -1;
  m_cnt_setLocalMatrix = 0;

  MPI_Comm_rank(*m_pcomm,&m_rank);

}

void Domain::SetMappingLoc2Glb(std::vector<int>& _l2g)
{
  m_l2g = _l2g;
  m_neq = (int) m_l2g.size();
  m_SetMappingGlb2Loc();
}

void Domain::m_SetMappingGlb2Loc()
{
  // maping DOF: global to local
  for (int i = 0; i < m_neq; i++)
    m_g2l[m_l2g[i]] = i;
}

void Domain::SetNeighboursRanks(const std::vector<int>& ranks) 
{

  m_neighboursRanks = ranks;
  
#if DBG > 0
  m_dbg_printNeighboursRanks();
#endif

}






void Domain::NumericAssemblingStiffnessAndRhs(
    std::vector<int>& glbIds,
    std::vector<double>& valLocK,
    std::vector<double>& valLocRHS)
{

  //TODO if patern does not change
  // use exixting triplets and replace Value() only

  if (m_cnt_setLocalMatrix==0)
  {
    m_rhs.resize(m_neq);
    m_rhs.setZero();
#if DBG > 0
    std::cout << "local numbering - start\n";
#endif
  }


  int neqLocal = glbIds.size();

  if (pow(neqLocal,2) != valLocK.size())
    std::runtime_error(__FILE__);

  int rowLocal(0);

  for (int row = 0; row < neqLocal; row++)
  {
    rowLocal = m_g2l[glbIds[row]];
    m_rhs(rowLocal) += valLocRHS[row];

    for (int col = 0; col < neqLocal; col++)
    {
      m_trK.push_back(
          T(rowLocal, m_g2l[glbIds[col]],
            valLocK[row + neqLocal * col]));
    }
#if DBG > 0
    std::cout << m_g2l[glbIds[row]] << ' ';
#endif
  }
#if DBG > 0
  std::cout << '\n';
#endif

  m_cnt_setLocalMatrix++;

}

void Domain::FinalizeStiffnessMatrixAndRhs()
{

  m_stiffnessMatrix.resize(m_neq, m_neq);
  m_stiffnessMatrix.setFromTriplets(m_trK.begin(), m_trK.end());
  m_stiffnessMatrix.makeCompressed();
  m_trK.clear();
  m_trK.shrink_to_fit();

#if DBG>1
  m_dbg_printStiffnessMatrix();
#endif
#if DBG>3
  m_dbg_printStiffnessMatrixSingularValues();
#endif

  //TODO add Neumann to m_rhs if exists

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

  int nInterf = m_neighboursRanks.size();
  if (nInterf == 0)
    std::runtime_error("undecomposed?");

  m_interfaces.resize(nInterf);
  for (int iR = 0; iR < nInterf; iR++)
    m_interfaces[iR].SetNeighbRank(m_neighboursRanks[iR]);


  
  //for (int iR = 0; iR < nInterf; iR++)
  for (auto& i_intfc : m_interfaces)
  {
    //int neighRank = m_interfaces[iR].GetNeighbRank();
    int neighRank = i_intfc.GetNeighbRank();
    int neq_neighb(0);


    hmpi.SendInt(&m_neq,1 , neighRank);
    hmpi.RecvInt(&neq_neighb,1 , neighRank);

    //m_interfaces[iR].SetNeighbNumbOfEqv(neq_neighb);
    i_intfc.SetNeighbNumbOfEqv(neq_neighb);
    
    std::vector<int> intersection(0);

    if (m_rank > neighRank)
    {
      // greather rank manages intersection
      std::vector<int> neighb_l2g(neq_neighb,0);
      hmpi.RecvInt(neighb_l2g.data(),neq_neighb , neighRank);
      intersection = tools::intersection(m_l2g, neighb_l2g);
    }
    else
    {
      hmpi.SendInt(m_l2g.data(),m_l2g.size(), neighRank);
    }

    int neqInterface(0);

    if (m_rank > neighRank)
    {
      neqInterface = intersection.size();
      hmpi.SendInt(&neqInterface,1, neighRank);
    }
    else
    {
      hmpi.RecvInt(&neqInterface,1, neighRank);
    }

    std::cout << "neqInterface= " << neqInterface << std::endl;




    if (m_rank > neighRank)
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
#if DBG > 3
    std::cout  << "#interf. dofs:  ";
    for (auto& ii : i_intfc.m_interfaceDOFs)
      std::cout << ii << ' ' ;
    std::cout<< std::endl;
#endif



  }
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


void Domain::m_dbg_printStiffnessMatrix()
{
  std::string fname = "matK" + std::to_string(m_rank) + ".txt";
  tools::printMatrix(m_stiffnessMatrix,fname);
}

void Domain::m_dbg_printStiffnessMatrixSingularValues()
{
  std::cout << "singular values of K: \n";
  auto singVals = linalg::svd0(m_stiffnessMatrix);
  for (int iv = 0; iv < singVals.size() ; iv++)
    std::cout << singVals[iv] << ' ';
  std::cout << '\n';
}


