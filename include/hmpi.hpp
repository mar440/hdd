#pragma once

#include <mpi.h>
#include <vector>



class Hmpi 
{

  public:
    Hmpi(MPI_Comm *_comm);
    ~Hmpi(){}

    int GetRank();
    int GetSize();


    MPI_Comm GetComm(){return m_comm;}
    void SendInt(void* buf,int count, int dest);
    void RecvInt(void* buf,int count, int dest);

    void IsendInt(void* buf,int count, int dest);
    void IrecvInt(void* buf,int count, int dest);


    void SendDbl(void* buf,int count, int dest);
    void RecvDbl(void* buf,int count, int dest);

    void IsendDbl(void* buf,int count, int dest);
    void IrecvDbl(void* buf,int count, int dest);

    void WaitRecv();
    void WaitSend();
    void Wait();

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

    void ScattervDbl(void* sendbuf, int* sendcnts, int* displs,
        void* recvbuf, int recvcnt, int root);

    void ScattervInt(void* sendbuf, int* sendcnts, int* displs,
        void* recvbuf, int recvcnt, int root);

    void BcastInt(void* buffer, int count, int root);

    void BcastDbl(void* buffer, int count, int root);

    void GlobalSum(double* in_out_buffer, int count);

    void GlobalInt(int* in_out_buffer, int count, MPI_Op operation);

    void GlobalInt(std::vector<int>& in_out_buffer, MPI_Op operation);

    void AlltoallInt(int* in_out_buffer,int in_buffer_size, int count);

    void ScatterInt(const void* sendbuf, int sendcount, void *recvbuf,
      int recvcount, int root);

  private:
    MPI_Comm m_comm;
    MPI_Status m_recv_status;
    MPI_Status m_send_status;
    MPI_Request m_send_request;
    MPI_Request m_recv_request;



};
