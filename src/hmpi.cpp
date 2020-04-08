#include "../include/hmpi.hpp"
#include <iostream>
#include <stdexcept>

#define TAG_SEND_INT 1000



Hmpi::Hmpi(MPI_Comm* _comm)
{
  m_pcomm = _comm;
}

void Hmpi::SendInt(void* buf,int count, int dest)
{
  MPI_Send(buf, count, MPI_INT, dest, TAG_SEND_INT, *m_pcomm );
}

void Hmpi::RecvInt(void* buf,int count, int dest)
{
  MPI_Status recv_status;
  MPI_Recv(buf, count,MPI_INT, dest ,TAG_SEND_INT,*m_pcomm,&recv_status);
}

void Hmpi::SendDbl(void* buf,int count, int dest)
{
  MPI_Send(buf, count, MPI_DOUBLE, dest, TAG_SEND_INT, *m_pcomm );
}

void Hmpi::RecvDbl(void* buf,int count, int dest)
{
  MPI_Status recv_status;
  MPI_Recv(buf, count,MPI_DOUBLE, dest ,TAG_SEND_INT,*m_pcomm,&recv_status);
}



void Hmpi::ReduceInt(const void *sendbuf, void *recvbuf, int count,
               MPI_Op op, int root)
{
  MPI_Reduce(sendbuf, recvbuf, count, MPI_INT, op, root, *m_pcomm);
}



void Hmpi::GatherInt(const void *sendbuf, int sendcount,
               void *recvbuf, int recvcount,int root)
{
  MPI_Gather(sendbuf, sendcount, MPI_INT, recvbuf, recvcount, MPI_INT,
               root, *m_pcomm);
}


void Hmpi::GathervInt(const void *sendbuf, int sendcount, 
                void *recvbuf, const int *recvcounts, const int *displs,
                int root)
{
  MPI_Gatherv(sendbuf, sendcount, MPI_INT, recvbuf, recvcounts,
                displs, MPI_INT, root, *m_pcomm);

}

void Hmpi::GathervDbl(const void *sendbuf, int sendcount, 
                void *recvbuf, const int *recvcounts, const int *displs,
                int root)
{
  MPI_Gatherv(sendbuf, sendcount, MPI_DOUBLE, recvbuf, recvcounts,
                displs, MPI_DOUBLE, root, *m_pcomm);

}

void Hmpi::BcastInt(void *buffer, int count, int root)
{
  MPI_Bcast(buffer,count, MPI_INT, root, *m_pcomm);
}

void Hmpi::BcastDbl(void *buffer, int count, int root)
{
  MPI_Bcast(buffer,count, MPI_DOUBLE, root, *m_pcomm);
}

void Hmpi::Barrier()
{
  MPI_Barrier(*m_pcomm);
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

