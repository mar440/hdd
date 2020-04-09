#pragma once

#include <mpi.h>

class Hmpi 
{

  public:
    Hmpi(MPI_Comm *_comm);
    ~Hmpi(){}

    void SendInt(void* buf,int count, int dest);
    void RecvInt(void* buf,int count, int dest);
    void SendDbl(void* buf,int count, int dest);
    void RecvDbl(void* buf,int count, int dest);
    void Sendrecv();
    void Barrier();

    void ReduceInt(const void *sendbuf, void *recvbuf, int count,
               MPI_Op op, int root);


    void GatherInt(const void *sendbuf, int sendcount,
               void *recvbuf, int recvcount,int root);

    void GathervInt(const void *sendbuf, int sendcount,
                void *recvbuf, const int *recvcounts, const int *displs,
                int root);

    void GathervDbl(const void *sendbuf, int sendcount,
                void *recvbuf, const int *recvcounts, const int *displs,
                int root);

    void ScattervDbl(void *sendbuf, int *sendcnts, int *displs,
        void *recvbuf, int recvcnt, int root);

    void BcastInt(void *buffer, int count, int root);

    void BcastDbl(void *buffer, int count, int root);

  private:
    MPI_Comm* m_pcomm;


};
