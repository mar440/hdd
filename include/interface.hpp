#pragma once
#include <vector>


class Interface
{
  public:
    Interface();
    ~Interface(){}

    int GetNeighbRank(){return m_neighb_rank;}
    void SetNeighbRank(int nr){m_neighb_rank = nr;}

    int GetNeighbNumbOfEqv(){return m_neq_neighb;}
    void SetNeighbNumbOfEqv(int _neq){m_neq_neighb = _neq;}

    int Get_neqInterface(){return m_neqInterface;}
//    void Set_neqInterface(int _neq){m_neqInterface= _neq;}
    
    std::vector<int> m_interfaceDOFs;

  private:
    int m_neighb_rank;
    int m_neq_neighb;
    int m_neqInterface;




};
