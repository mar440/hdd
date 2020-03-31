#pragma once

#include <mpi.h>

class Hmpi 
{

  public:
    Hmpi(MPI_Comm *_comm);
    ~Hmpi(){}

    void SendInt(void* buf,int count, int dest);
    void RecvInt(void* buf,int count, int dest);
    void Barrier();

  private:
    MPI_Comm* m_pcomm;


};
