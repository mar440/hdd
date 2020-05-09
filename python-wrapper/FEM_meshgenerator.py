from __future__ import division
from __future__ import print_function
import numpy as np


def MeshGenerator(example_options):

    Nx = example_options['Nx']
    Ny = example_options['Ny']
    Nz = example_options['Nz']

    nx = example_options['nx'] * example_options['Nx']
    ny = example_options['ny'] * example_options['Ny']
    nz = example_options['nz'] * example_options['Nz']

    Lx = example_options['Lx']
    Ly = example_options['Ly']
    Lz = example_options['Lz']

    nnods = (nx+1)*(ny+1)*(nz+1)
    nnods_h = "nx2 array - nodes coordinates"
    nelem =  nx*ny*nz
    coordinates = np.zeros((nnods,3),np.float64)
    elements    = np.zeros((nelem,8),np.int32)   
    edges       = np.zeros((2*nx*ny+2*nx*nz+2*ny*nz,7),np.int32) 
    # coordinates
    Z,Y,X = np.mgrid[0:nz+1,0:ny+1,0:nx+1] 
    coordinates[:,0] = X.reshape(nnods)*(float(Lx)/nx)
    coordinates[:,1] = Y.reshape(nnods)*(float(Ly)/ny)
    coordinates[:,2] = Z.reshape(nnods)*(float(Lz)/nz)
    # elements 
    inodes  = np.arange(0,(nx+1)*(ny+1)*(nz+1),1)
    inodes  = inodes.reshape(nz+1,ny+1,nx+1)
    elements[0:(nx*ny*nz),0]   = np.reshape(inodes[0:-1,    0:-1,   0:-1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),1]   = np.reshape(inodes[0:-1,    0:-1,   1:nx+1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),2]   = np.reshape(inodes[0:-1,    1:ny+1, 1:nx+1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),3]   = np.reshape(inodes[0:-1,    1:ny+1, 0:-1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),4]   = np.reshape(inodes[1:nz+1,  0:-1,   0:-1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),5]   = np.reshape(inodes[1:nz+1,  0:-1,   1:nx+1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),6]   = np.reshape(inodes[1:nz+1,  1:ny+1, 1:nx+1],(nz)*(ny)*(nx),1) 
    elements[0:(nx*ny*nz),7]   = np.reshape(inodes[1:nz+1,  1:ny+1, 0:-1],(nz)*(ny)*(nx),1)  
    # edges
    cnt = 0
    # E0-------------------------------------------
    E0  = inodes[0,:,:]
    E00 = np.reshape(E0[0:ny,0:nx],nx*ny)
    E01 = np.reshape(E0[0:ny,1:nx+1],nx*ny)
    E02 = np.reshape(E0[1:ny+1,1:nx+1],nx*ny)
    E03 = np.reshape(E0[1:ny+1,0:nx],nx*ny)
    ind=np.arange(0,nx*ny)+cnt
    cnt = cnt + nx*ny
    edges[ind,0] = E00
    edges[ind,1] = E01
    edges[ind,2] = E02
    edges[ind,3] = E03
    edges[ind,4] = 0
    # E1-------------------------------------------
    E1  = inodes[-1,:,:] 
    E10 = np.reshape(E1[0:ny,0:nx],nx*ny)
    E11 = np.reshape(E1[0:ny,1:nx+1],nx*ny)
    E12 = np.reshape(E1[1:ny+1,1:nx+1],nx*ny)
    E13 = np.reshape(E1[1:ny+1,0:nx],nx*ny) 
    ind=np.arange(0,(nx*ny))+cnt
    cnt = cnt + nx*ny
    edges[ind,0] = E10
    edges[ind,1] = E13
    edges[ind,2] = E12
    edges[ind,3] = E11
    edges[ind,4] = 1
    # E2-------------------------------------------
    E2  = inodes[:,0,:] 
    E20 = np.reshape(E2[1:nz+1,0:nx],nx*nz)
    E21 = np.reshape(E2[0:nz,0:nx],nx*nz)
    E22 = np.reshape(E2[0:nz,1:nx+1],nx*nz)
    E23 = np.reshape(E2[1:nz+1,1:nx+1],nx*nz) 
    ind=np.arange(0,nx*nz)+cnt
    cnt = cnt + nx*nz
    edges[ind,0] = E20
    edges[ind,1] = E23
    edges[ind,2] = E22
    edges[ind,3] = E21       
    edges[ind,4] = 2     
    # E3-------------------------------------------
    E3  = inodes[:,-1,:] 
    E30 = np.reshape(E3[0:nz,0:nx],nx*nz)
    E31 = np.reshape(E3[1:nz+1,0:nx],nx*nz)
    E32 = np.reshape(E3[1:nz+1,1:nx+1],nx*nz)
    E33 = np.reshape(E3[0:nz,1:nx+1],nx*nz) 
    ind=np.arange(0,nx*nz)+cnt
    cnt = cnt + nx*nz
    edges[ind,0] = E30
    edges[ind,1] = E33
    edges[ind,2] = E32
    edges[ind,3] = E31
    edges[ind,4] = 3
    # E4-------------------------------------------
    E4  = inodes[:,:,0] 
    E40 = np.reshape(E4[0:nz,0:ny],ny*nz)
    E41 = np.reshape(E4[0:nz,1:ny+1],ny*nz)
    E42 = np.reshape(E4[1:nz+1,1:ny+1],ny*nz)
    E43 = np.reshape(E4[1:nz+1,0:ny],ny*nz) 
    ind=np.arange(0,ny*nz)+cnt
    cnt = cnt + ny*nz
    edges[ind,0] = E40
    edges[ind,1] = E41
    edges[ind,2] = E42
    edges[ind,3] = E43         
    edges[ind,4] = 4   
    # E5-------------------------------------------
    E5  = inodes[:,:,-1] 
    E50 = np.reshape(E5[1:nz+1,0:ny],ny*nz)
    E51 = np.reshape(E5[1:nz+1,1:ny+1],ny*nz)
    E52 = np.reshape(E5[0:nz,1:ny+1],ny*nz)
    E53 = np.reshape(E5[0:nz,0:ny],ny*nz) 
    ind=np.arange(0,ny*nz)+cnt
    cnt = cnt + ny*nz
    edges[ind,0] = E50
    edges[ind,1] = E51
    edges[ind,2] = E52
    edges[ind,3] = E53
    edges[ind,4] = 5             

    if Nx*Ny*Nz>1:
        #a[ix_([1,3,4],[0,2])]
        X=np.zeros((nx,ny,nz),np.int32)
        c=0
        for i in range(Nx):
            for j in range(Ny):
                for k in range(Nz):
                    nNx = int(nx/Nx) 
                    nNy = int(ny/Ny)
                    nNz = int(nz/Nz)
                    in1 = np.arange((i)*nNx,(i+1)*nNx)
                    in2 = np.arange((j)*nNy,(j+1)*nNy)
                    in3 = np.arange((k)*nNz,(k+1)*nNz)
                    X[np.ix_(in1,in2,in3)]=c
                    c=c+1

        subs0 = X.reshape(nx*ny*nz,order='F')
    else:
        subs0 = np.zeros(1)

    return {'coordinates':coordinates,'elements': elements, 'edges':edges,'decomp':subs0}

    #################################################################
    #                                                               #
    #                         A z                                   #  
    #                         |                                     #    
    #                         |            L4                       #
    #                         |_ _ _ _ _ _ _                        #
    #                        /     L1      /|                       #  
    #                       /_ _ _ _ _ _  / |                       #
    #                      |      |      |  |                       #  
    #                    L2|      |      |L3|                       #   
    #                      |_ _ _ |_ _ _ |  |                       #  
    #                      |    L5|      |  |-------> y             #
    #                      |      |      | /                        #
    #                      |_ _ _ |_ _ _ |/                         #
    #                     /                                         #  
    #                    /       L0                                 #
    #                   /                                           #  
    #                  v  x                                         #
    #                                                               #
    #################################################################   



