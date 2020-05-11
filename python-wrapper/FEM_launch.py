from __future__ import division
from __future__ import print_function


import numpy as np
from pythonHdd import PyHddApi
from mpi4py import MPI

import sys
import trace

import FEM_meshgenerator as m3d
import FEM_assembler as assmbl

mpicomm = MPI.COMM_WORLD
mpirank = mpicomm.rank

launchWithTraces = False

def runParallel():
    nx = 4
    ny = 4
    nz = 4

    Nx = 3
    Ny = 3
    Nz = 3

    Lx = Ly = Lz = 1.0


    if mpirank==0:print("start ... Python")
    hdd= PyHddApi(mpicomm)
    hdd.Message("Beginning ...")
    hdd.ParseJsonFile("hddConf.json")

    mesh = m3d.MeshGenerator({'Nx':Nx, 'Ny':Ny , 'Nz':Nz,
                'nx':nx , 'ny':ny , 'nz' : nz,
                'Lx':Lx , 'Ly':Ly , 'Lz' : Lz })

    dimension = 3
    Ex = 1.1
    mu = 0.3

    coordinates = mesh['coordinates']
    elements    = mesh['elements']
    #edges       = mesh['edges']
    decomp      = mesh['decomp']
    nElements   = elements.shape[0]
    volF        = np.ones(3)

    DirPts = np.where(mesh['coordinates'][:,2] == 0)[0]
    nDP = DirPts.shape[0]

    dirichletDofs = np.zeros(3 * nDP,dtype=np.int32)
    dirichletDofs[0::3] = 3 * DirPts + 0
    dirichletDofs[1::3] = 3 * DirPts + 1
    dirichletDofs[2::3] = 3 * DirPts + 2


    for i in range(nElements):
        if (decomp[i] == mpirank):
            ielem = elements[i,:]
            nP =  ielem.shape[0]
            neqLoc = nP * dimension
            indx = np.zeros(neqLoc,dtype=np.int32)
            indx[0*nP:1*nP:] = ielem * 3 + 0
            indx[1*nP:2*nP:] = ielem * 3 + 1
            indx[2*nP:3*nP:] = ielem * 3 + 2
            hdd.SymbolicAssembling(indx)
    ###
    hdd.SetDirichletDOFs(dirichletDofs)
    hdd.FinalizeSymbolicAssembling()
    #
    for i in range(nElements):
        if (decomp[i] == mpirank):
            ielem = elements[i,:]
            coordinatesLoc = coordinates[ielem,]
            nP =  ielem.shape[0]
            neqLoc = nP * dimension
            indx = np.zeros(neqLoc,dtype=np.int32)
            indx[0*nP:1*nP:] = ielem * 3 + 0
            indx[1*nP:2*nP:] = ielem * 3 + 1
            indx[2*nP:3*nP:] = ielem * 3 + 2

            localElemData = assmbl.stima3e8({'coordinatesLoc':coordinatesLoc,
                'Ex':Ex, 'mu':mu, 'volF':volF})
            Ke = localElemData['Ke'].copy()
            Fe = localElemData['fe'].copy()
            hdd.NumericAssembling(indx,Ke,Fe)

    hdd.FinalizeNumericAssembling()


    solution = hdd.Solve()

#    print(solution)

    if mpirank == 0: print("end ... Python")



if launchWithTraces:
    # define Trace object: trace line numbers at runtime, exclude some modules
    tracer = trace.Trace(
        ignoredirs=[sys.prefix, sys.exec_prefix],
        ignoremods=[
            'inspect', 'contextlib', '_bootstrap',
            '_weakrefset', 'abc', 'posixpath', 'genericpath', 'textwrap'
        ],
        trace=1,
        count=0)

    # by default trace goes to stdout
    # redirect to a different file for each processes
    sys.stdout = open('trace_{:04d}.txt'.format(mpirank), 'w')

    tracer.runfunc(runParallel)

else:
    runParallel()
