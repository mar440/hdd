#include "../include/hmpi.hpp"
#include "../include/types.hpp"
#include "../include/hddTime.hpp"

#include <iostream>
#include <stdexcept>
#include <vector>

#define TAG_SEND_INT 1000
#define TAG_SEND_DBL 2000



Hmpi::Hmpi(MPI_Comm* _pcomm)
{
  m_comm = *_pcomm;
}

int Hmpi::GetRank()
{
  int _rank;
  MPI_Comm_rank(m_comm,&_rank);
  return _rank;
}

int Hmpi::GetSize()
{
  int _size;
  MPI_Comm_size(m_comm,&_size);
  return _size;
}

void Hmpi::SendInt(void* buf,int count, int dest)
{
  MPI_Send(buf, count, MPI_INT, dest, TAG_SEND_INT, m_comm );
}

void Hmpi::IsendInt(void* buf,int count, int dest)
{
  MPI_Isend(buf, count, MPI_INT, dest, TAG_SEND_INT, m_comm, &m_send_request);
}

void Hmpi::RecvInt(void* buf,int count, int dest)
{
  MPI_Status recv_status;
  MPI_Recv(buf, count,MPI_INT, dest ,TAG_SEND_INT,m_comm,&recv_status);
}

void Hmpi::IrecvInt(void* buf,int count, int dest)
{
  MPI_Irecv(buf,count,MPI_INT, dest, TAG_SEND_INT, m_comm,&m_recv_request);
}

void Hmpi::SendDbl(void* buf,int count, int dest)
{
  MPI_Send(buf, count, MPI_DOUBLE, dest, TAG_SEND_DBL, m_comm );
}

void Hmpi::IsendDbl(void* buf,int count, int dest)
{
  MPI_Isend(buf, count, MPI_DOUBLE, dest, TAG_SEND_DBL, m_comm, &m_send_request );
}

void Hmpi::RecvDbl(void* buf,int count, int dest)
{
  MPI_Status recv_status;
  MPI_Recv(buf, count,MPI_DOUBLE, dest ,TAG_SEND_DBL,m_comm,&recv_status);
}

void Hmpi::IrecvDbl(void* buf,int count, int dest)
{
  MPI_Irecv(buf,count,MPI_DOUBLE, dest, TAG_SEND_DBL, m_comm,&m_recv_request);
}

void Hmpi::Wait()
{
  WaitRecv();
  WaitSend();
}

void Hmpi::WaitRecv()
{
  MPI_Wait(&m_recv_request, &m_recv_status);
}

void Hmpi::WaitSend()
{
  MPI_Wait(&m_send_request, &m_send_status);
}


void Hmpi::ReduceInt(const void *sendbuf, void *recvbuf, int count,
               MPI_Op op, int root)
{
  MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, op, root, m_comm);
}

void Hmpi::GatherInt(const void *sendbuf, int sendcount,
               void *recvbuf, int recvcount,int root)
{
  MPI_Gather(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT,
               root, m_comm);
}

void Hmpi::GathervInt(const void *sendbuf, int sendcount, 
                void *recvbuf, const int *recvcounts, const int *displs,
                int root)
{
  MPI_Gatherv(sendbuf, sendcount, MPI_INT, recvbuf, recvcounts,
                displs, MPI_INT, root, m_comm);
}

void Hmpi::GathervDbl(const void *sendbuf, int sendcount, 
                void *recvbuf, const int *recvcounts, const int *displs,
                int root)
{
  MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcounts,
                displs, MPI_DOUBLE, root, m_comm);
}

void Hmpi::ScattervDbl(void *sendbuf, int *sendcnts, int *displs,
  void *recvbuf, int recvcnt, int root)
{

  MPI_Scatterv(sendbuf,sendcnts,displs,MPI_DOUBLE,
      recvbuf,recvcnt, MPI_DOUBLE,root,m_comm);
}

void Hmpi::ScattervInt(void *sendbuf, int *sendcnts, int *displs,
  void *recvbuf, int recvcnt, int root)
{

  MPI_Scatterv(sendbuf,sendcnts,displs,MPI_INT,
      recvbuf,recvcnt, MPI_INT,root,m_comm);
}

void Hmpi::BcastInt(void *buffer, int count, int root)
{
  MPI_Bcast(buffer,count, MPI_INT, root, m_comm);
}

void Hmpi::BcastDbl(void *buffer, int count, int root)
{
  MPI_Bcast(buffer,count, MPI_DOUBLE, root, m_comm);
}


void Hmpi::GlobalSum(double * in_out_buffer, int count)
{
  std::vector<double> send_buffer(in_out_buffer,  in_out_buffer + count);

  MPI_Allreduce(send_buffer.data(), in_out_buffer, count, MPI_DOUBLE,
      MPI_SUM,m_comm);
}

void Hmpi::GlobalInt(int* in_out_buffer, int count, MPI_Op operation)
{
  std::vector<int> send_buffer(in_out_buffer,  in_out_buffer + count);

  MPI_Allreduce(send_buffer.data(), in_out_buffer, count, MPI_INT,
      operation, m_comm);
}

void Hmpi::GlobalInt(std::vector<int>& in_out_buffer, MPI_Op operation)
{

  HDDTRACES
  std::vector<int> tmp(in_out_buffer);
  HDDTRACES
  int mpimsg = MPI_Allreduce(tmp.data(), in_out_buffer.data(), in_out_buffer.size(), 
      MPI_INT, operation, m_comm);
  std::cout << "mpimsg: " << mpimsg << std::endl;
  HDDTRACES
}

void Hmpi::Barrier()
{
  MPI_Barrier(m_comm);
}


void Hmpi::AlltoallInt(int* in_out_buffer, int in_buffer_size, int count)
{

  std::vector<int> send_buffer(in_out_buffer,  in_out_buffer + in_buffer_size);

  MPI_Alltoall(send_buffer.data(),count,MPI_INT,in_out_buffer,count,MPI_INT,
      m_comm);

}

void Hmpi::ScatterInt(const void* sendbuf, int sendcount, void *recvbuf, 
    int recvcount, int root)
{
    MPI_Scatter(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT,
        root, m_comm);
}



void Hmpi::Sendrecv()
{
  std::runtime_error("Sendrecv is not implemented");
//
//
//  int MPI_Sendrecv(const void *sendbuf, int sendcount, MPI_Datatype sendtype,
//                int dest, int sendtag,
//                void *recvbuf, int recvcount, MPI_Datatype recvtype,
//                int source, int recvtag,
//                MPI_Comm comm, MPI_Status *status)

}

