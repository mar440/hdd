#include "../include/hmpi.hpp"

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


void Hmpi::Barrier()
{
  MPI_Barrier(*m_pcomm);
}

////     MPI_Isend(void* buf, int count, MPI_Datatype datatype, 
////	     int dest, int tag, MPI_Comm comm, 
////	     MPI_Request *request);
//      int sizeOfMyBuffer = m_neq;
//      int sizeOfNeighBuffer(0);
//
//      std::cout << "sizeOfNeighBuffer: " << sizeOfNeighBuffer << std::endl;
//
////      MPI_Request send_request, recv_request;
//
////      MPI_Isend(&sizeOfMyBuffer, 1, MPI_INT, neighRank, tag, *m_pcomm, &send_request);
//
////      MPI_Irecv (&sizeOfNeighBuffer, 1, MPI_INT, neighRank, tag, *m_pcomm, &recv_request);
//
//      MPI_Send(&sizeOfMyBuffer, 1, MPI_INT, neighRank, tag, *m_pcomm );
////      MPI_Send(&sizeOfMyBuffer, 1, MPI_INT, neighRank, tag, *m_pcomm);
//      MPI_Recv(&sizeOfNeighBuffer,1,MPI_INT,neighRank,tag,*m_pcomm,&recv_status);
//
//      std::cout << "sizeOfNeighBuffer: " << sizeOfNeighBuffer << std::endl;

