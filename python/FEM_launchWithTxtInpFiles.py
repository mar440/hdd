from __future__ import division
from __future__ import print_function

import numpy as np
from pythonHdd import PyHddApi
from mpi4py import MPI

comm = MPI.COMM_WORLD
mpirank = comm.rank
mpisize = comm.size

if mpirank==0: print("start ... Python")

wdir = "/home/mar440/WorkSpace/test/"
fileBuffer = [] # to 'close' all files in the end

fnameGlbDirichlet   = "dmpGlbDirichletIds_"+str(mpirank)+".txt"
fGlbDirichlet = open(wdir+fnameGlbDirichlet)
fileBuffer.append(fGlbDirichlet)

fnameGlbDOFs        = "dmpGlbDOFs_"+str(mpirank)+".txt"
fGlbDOFs = open(wdir+fnameGlbDOFs)
fileBuffer.append(fGlbDOFs)

fnameLocLinOper     = "dmpLocLinOperator_"+str(mpirank)+".txt"
fLocLinOper = open(wdir+fnameLocLinOper)
fileBuffer.append(fLocLinOper)

fnameLocRHS         = "dmpLocRHS_"+str(mpirank)+".txt"
fLocRHS= open(wdir+fnameLocRHS);
fileBuffer.append(fLocRHS)


hdd = PyHddApi(comm)
hdd.Message("Beginning ...")
hdd.ParseJsonFile(wdir+"hddConf.json")
GlbDOFs_readlines = fGlbDOFs.readlines()

nElements = len(GlbDOFs_readlines)
for iline in GlbDOFs_readlines:
    GDOFS = np.array(iline.split(),dtype=np.int32)
    hdd.SymbolicAssembling(GDOFS)
DirDOFs = np.array(fGlbDirichlet.readline().split(),dtype=np.int32)
hdd.SetDirichletDOFs(DirDOFs)
hdd.FinalizeSymbolicAssembling()

LocLinOper_readlines = fLocLinOper.readlines()
LocRHS_readlines = fLocRHS.readlines()

for cnt in range(nElements):
    GDOFS = \
        np.array(LocLinOper_readlines[2*cnt].split(),
                dtype=np.int32)
    valsK = \
        np.array(LocLinOper_readlines[2*cnt +1].split(),
                dtype=np.float64)
    valsF = \
        np.array(LocRHS_readlines[2*cnt +1].split(),
                dtype=np.float64)
    neqLoc = GDOFS.shape[0]
    valsK = valsK.reshape((neqLoc,neqLoc))
    hdd.NumericAssembling(GDOFS,valsK,valsF)

hdd.FinalizeNumericAssembling()
solution = hdd.Solve()

np.savetxt(wdir+"sol_"+str(mpirank)+".txt",solution)
if mpirank == 0: print("end ... Python")

for ff in fileBuffer:
    ff.close()

