#pragma once
#include <vector>

class Domain;
class InterfaceOperatorB;

class Interface
{
  friend Domain;
  friend InterfaceOperatorB;
  public:
    Interface();
    ~Interface(){}

    int GetNeighbRank(){return m_neighb_rank;}
    void SetNeighbRank(int nr){m_neighb_rank = nr;}

    int GetNeighbNumbOfEqv(){return m_neq_neighb;}
    void SetNeighbNumbOfEqv(int _neq){m_neq_neighb = _neq;}

    int GetNeighbDefect(){return m_neighb_defect;}
    void SetNeighbDefect(int _defect){m_neighb_defect= _defect;}


    int Get_neqInterface(){return m_neqOnInterface;}

  protected:
    std::vector<int> m_interfaceDOFs;

  private:
    int m_neighb_rank;
    int m_neq_neighb;
    int m_neqOnInterface;
    int m_neighb_defect;
    int m_offset;
};
