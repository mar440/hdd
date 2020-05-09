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
fileBuffer = [] # for to easy 'close' all txt-files in the end

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
hdd.Message("Beginning ...") # pass any message to hdd-output file
hdd.ParseJsonFile(wdir+"hddConf.json")
GlbDOFs_readlines = fGlbDOFs.readlines()

nElements = len(GlbDOFs_readlines)
for iline in GlbDOFs_readlines:
    eqIndx = np.array(iline.split(),dtype=np.int32)
    hdd.SymbolicAssembling(eqIndx)
DirDOFs = np.array(fGlbDirichlet.readline().split(),dtype=np.int32)
hdd.SetDirichletDOFs(DirDOFs)
hdd.FinalizeSymbolicAssembling()

LocLinOper_readlines = fLocLinOper.readlines()
LocRHS_readlines = fLocRHS.readlines()

for cnt in range(nElements):
    eqIndx = \
        np.array(LocLinOper_readlines[2*cnt].split(),
                dtype=np.int32)
    locK = \
        np.array(LocLinOper_readlines[2*cnt +1].split(),
                dtype=np.float64)
    locRHS = \
        np.array(LocRHS_readlines[2*cnt +1].split(),
                dtype=np.float64)
    neqLoc = eqIndx.shape[0]
    locK = locK.reshape((neqLoc,neqLoc))
    hdd.NumericAssembling(eqIndx,locK,locRHS)

hdd.FinalizeNumericAssembling()
solution = hdd.Solve()

np.savetxt(wdir+"sol_"+str(mpirank)+".txt",solution, fmt='%.18e')
if mpirank == 0: print("end ... Python")

for ff in fileBuffer:
    ff.close()

