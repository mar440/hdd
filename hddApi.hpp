#pragma once

#include <string>
#include "include/data.hpp"
#include <mpi.h>



struct HddApi
{
  public:
    HddApi(MPI_Comm* pcomm) : m_data(pcomm){}
    HddApi(MPI_Comm comm) : m_data(&comm){}
    ~HddApi(){}
    void ParseJsonFile(std::string);
    void Test(std::string str){std::cout<<str  << " in c++\n";}
    void SymbolicAssembling(int vals[], int size);
    void SetDirichletDOFs(int vals[], int size);
    int FinalizeSymbolicAssembling();
    void NumericAssembling(int glbIds[],
        double valK[], double valF[], int nDOF);
    void FinalizeNumericAssembling();
    void Solve(double solution[], int nDofOnSubdomain, int ncols = 1);
    void Finalize();

  private:
    Data m_data;

};
