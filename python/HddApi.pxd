# distutils: language = c++

from libcpp.string cimport string
from mpi4py cimport libmpi

cdef extern from "src/hddApi.cpp":
    pass

# Declare the class with cdef
cdef extern from "hddApi.hpp":
    cdef cppclass HddApi:
        HddApi(libmpi.MPI_Comm _comm) except +
        void ParseJsonFile(string)
        void SymbolicAssembling(int vals[], int size)
        void SetDirichletDOFs(int vals[], int size)
        void Test(string)
        int FinalizeSymbolicAssembling()
        void NumericAssembling(int glbIds[],
            double valK[], double valF[], int nDOF)
        void FinalizeNumericAssembling()
        void Solve(double solution[], int nDofOnSubdomain, int ncols )
        void Finalize()
