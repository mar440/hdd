from __future__ import division
from __future__ import print_function


import numpy as np
from pythonHdd import PyHddApi
from mpi4py import MPI

import FEM_meshgenerator as m3d
import FEM_assembler as assmbl


nx = ny = nz = 4

Nx = 1
Ny = 2
Nz = 2

Lx = 1.0
Ly = 1.0
Lz = 1.0

comm = MPI.COMM_WORLD
rank = comm.rank

if True:

    print("start ... Python")
    hdd= PyHddApi(comm)
    hdd.Message("Beginning ...")
if True:
    hdd.ParseJsonFile("/home/mar440/WorkSpace/hfls/mpfeti_separate/python/hddConf.json")

    mesh = m3d.MeshGenerator({'Nx':Nx, 'Ny':Ny , 'Nz':Nz,
                'nx':nx , 'ny':ny , 'nz' : nz,
                'Lx':Lx , 'Ly':Ly , 'Lz' : Lz })


    dimension = 3; Ex = 1.1; mu = 0.3

    coordinates = mesh['coordinates']
    elements    = mesh['elements']
    edges       = mesh['edges']
    decomp      = mesh['decomp']
    nElements = elements.shape[0]
    volF        = np.ones(3)

    DirichletFace = 0
    DirPts = \
        np.unique(edges[edges[:,4]==DirichletFace,0:4].ravel())
    nDP = DirPts.shape[0]

    dirichletDofs = np.zeros(3 * nDP,dtype=np.int32)
    dirichletDofs[0::3] = 3 * DirPts + 0
    dirichletDofs[1::3] = 3 * DirPts + 1
    dirichletDofs[2::3] = 3 * DirPts + 2


    for i in range(nElements):
        if (decomp[i] == rank):
            ielem = elements[i,:]
            nP =  ielem.shape[0]
            neqLoc = nP * dimension
            # ul = [ul1x,ul2x,ul3x,...., ul1y,ul2y,ul3y,.....]
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
        if (decomp[i] == rank):
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
            Ke = localElemData['Ke']
            Fe = localElemData['fe']
            hdd.NumericAssembling(indx,Ke,Fe)

    hdd.FinalizeNumericAssembling()


    solution = hdd.Solve()

    print(solution)

    if rank == 0: print("end ... Python")
