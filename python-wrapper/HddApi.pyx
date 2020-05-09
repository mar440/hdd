# distutils: language = c++

from HddApi cimport HddApi as cHddApi
from libcpp.string cimport string

cimport mpi4py.MPI as MPI
import numpy as np
cimport numpy as np


cdef class PyHddApi:
    cdef cHddApi *thisptr
    cdef int __neqOnSubdomain__
    def __cinit__(self,MPI.Comm comm):
        self.thisptr = new cHddApi(comm.ob_mpi)
    def __dealloc__(self):
        self.thisptr.Finalize()
        del self.thisptr
    def Message(self, str):
        self.thisptr.Test(str)
    def ParseJsonFile(self,str):
        self.thisptr.ParseJsonFile(str)
    def SymbolicAssembling(self,np.ndarray[int, ndim=1] indxs not None):
        cdef int n
        n = indxs.shape[0]
        self.thisptr.SymbolicAssembling(&indxs[0],n)
    def SetDirichletDOFs(self,np.ndarray[int, ndim=1,mode="c"] indxs not None):
        cdef int n
        n = indxs.shape[0]
        self.thisptr.SetDirichletDOFs(&indxs[0],n)
    def FinalizeSymbolicAssembling(self):
        cdef int neqSub
        neqSub = self.thisptr.FinalizeSymbolicAssembling()
        self.__neqOnSubdomain__ = neqSub
    def NumericAssembling(self, np.ndarray[int, ndim=1,mode="c"] glbIds not None,\
                            np.ndarray[double, ndim=2,mode="c"] valK not None,\
                            np.ndarray[double, ndim=1,mode="c"] valF not None):
        cdef nI
        cdef nK
        cdef mK
        cdef nF
        nI = glbIds.shape[0]
        nK,mK = valK.shape[0],valK.shape[1]
        nF = valF.shape[0]
        self.thisptr.NumericAssembling(&glbIds[0], &valK[0,0],&valF[0], nI)
    def FinalizeNumericAssembling(self):
        self.thisptr.FinalizeNumericAssembling()
    def GetNeqSub(self):
        return self.__neqOnSubdomain__

    def Solve(self):
        cdef np.ndarray[double, ndim=1, mode="c"] sol = np.zeros(self.__neqOnSubdomain__)
        cdef double* psol = &sol[0]
        cdef int ncols
        ncols = 1
        self.thisptr.Solve(psol, self.__neqOnSubdomain__,ncols)
        return sol



