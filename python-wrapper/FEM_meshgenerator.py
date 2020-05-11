from __future__ import division
from __future__ import print_function
import numpy as np



def MeshGenerator(example_options):

    Nx = example_options['Nx']
    Ny = example_options['Ny']
    Nz = example_options['Nz']

    nxd = example_options['nx']# * example_options['Nx']
    nyd = example_options['ny']# * example_options['Ny']
    nzd = example_options['nz']# * example_options['Nz']

    Lx = example_options['Lx']
    Ly = example_options['Ly']
    Lz = example_options['Lz']

    x0 = y0 = z0 = 0.0

#def create_fem3d(   Lx = 1., Ly = 1., Lz = 1., x0 = 0.0, y0 = 0.0, z0 = 0.0,\
#                    nxd = 3, nyd = 3, nzd = 3,  Nx = 2, Ny = 2, Nz = 2):

#    Lx = 1.; Ly = 1.; Lz = 1.; nxd = 3;nyd = 3; nzd = 3; Nx = 2; Ny = 2; Nz = 2
#    x0 = 0.0; y0 = 0.0; z0 = 0.0;

    if Lx == 0 and Ly == 0 and Lz == 0:
        nEl = 0
        nNod = 0
        elements=np.zeros((nEl,8),dtype=np.int32)
        coordinates=np.zeros((nNod,3),dtype=np.float64)
        MaterialId = np.ones(nEl,dtype = np.int32)
        MaterialId = np.ones(nEl,dtype = np.int32)
        PartitionId = np.zeros(nEl,dtype = np.int32)
        FormulationId = np.ones(nEl,dtype = np.int32) 
    else:
#        return elements, coordinates, coarse_element_id, MaterialId, PartitionId, FormulationId
     
        nx = nxd * Nx
        ny = nyd * Ny
        nz = nzd * Nz
        nnods = (nx+1)*(ny+1)*(nz+1)
        nelem =  nx*ny*nz
        coordinates = np.zeros((nnods,3),dtype = np.float64)
        elements    = np.zeros((nelem,8),dtype = np.int32)    
        
        
        cnt = 0
        nxy = (nx + 1)*(ny + 1)
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    tmp = np.array([i + 0 + (j + 0) * (nx + 1),
                                    i + 1 + (j + 0) * (nx + 1),
                                    i + 1 + (j + 1) * (nx + 1),
                                    i + 0 + (j + 1) * (nx + 1)])
                    elements[cnt,:] = np.concatenate((tmp + k * nxy, tmp + (k + 1) * nxy))
                    cnt += 1
        
        dx = Lx/nx
        dy = Ly/ny
        dz = Lz/nz
        cnt = 0
        for k in range(nz+1):
            for j in range(ny+1):
                for i in range(nx+1):
                    coordinates[cnt,:] = np.array([i*dx + x0 , j*dy + y0 , k*dz + z0])
                    cnt += 1
        if Nx*Ny*Nz > 1:
            X=np.zeros((nx,ny,nz), dtype = np.int32)
            c=0
            for k in range(Nz):
                for j in range(Ny):
                    for i in range(Nx):
                        in1 = np.arange((i)*nxd,(i+1)*nxd)
                        in2 = np.arange((j)*nyd,(j+1)*nyd)
                        in3 = np.arange((k)*nzd,(k+1)*nzd)
                        X[np.ix_(in1,in2,in3)]=c
                        c+=1 
            PartitionId = X.reshape(nelem, order='F') 
        else:
            PartitionId = np.zeros(nelem, dtype = np.int32) 
 
        MaterialId = np.ones(nelem, dtype = np.int32)
        FormulationId = 30 * np.ones(nelem, dtype = np.int32 )
         
        
    #    elements = 0
    #    coordinates = 0 
    #    PartitionId = 0
    #    MaterialId = 0
    #    FormulationId = 0

    return {'coordinates':coordinates,'elements': elements,'decomp':PartitionId}
#    return elements, coordinates, PartitionId, MaterialId, FormulationId 


#def MeshGenerator(example_options):
#
#    Nx = example_options['Nx']
#    Ny = example_options['Ny']
#    Nz = example_options['Nz']
#
#    nx = example_options['nx'] * example_options['Nx']
#    ny = example_options['ny'] * example_options['Ny']
#    nz = example_options['nz'] * example_options['Nz']
#
#    Lx = example_options['Lx']
#    Ly = example_options['Ly']
#    Lz = example_options['Lz']
#
#    nnods = (nx+1)*(ny+1)*(nz+1)
#    nnods_h = "nx2 array - nodes coordinates"
#    nelem =  nx*ny*nz
#    coordinates = np.zeros((nnods,3),np.float64)
#    elements    = np.zeros((nelem,8),np.int32)   
##    edges       = np.zeros((2*nx*ny+2*nx*nz+2*ny*nz,7),np.int32) 
#    edges = []
#    # coordinates
#    Z,Y,X = np.mgrid[0:nz+1,0:ny+1,0:nx+1] 
#    coordinates[:,0] = X.reshape(nnods)*(np.float64(Lx)/nx)
#    coordinates[:,1] = Y.reshape(nnods)*(np.float64(Ly)/ny)
#    coordinates[:,2] = Z.reshape(nnods)*(np.float64(Lz)/nz)
#    # elements 
#    inodes  = np.arange(0,(nx+1)*(ny+1)*(nz+1),1)
#    inodes  = inodes.reshape(nz+1,ny+1,nx+1)
#    elements[0:(nx*ny*nz),0]   = np.reshape(inodes[0:-1,    0:-1,   0:-1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),1]   = np.reshape(inodes[0:-1,    0:-1,   1:nx+1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),2]   = np.reshape(inodes[0:-1,    1:ny+1, 1:nx+1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),3]   = np.reshape(inodes[0:-1,    1:ny+1, 0:-1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),4]   = np.reshape(inodes[1:nz+1,  0:-1,   0:-1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),5]   = np.reshape(inodes[1:nz+1,  0:-1,   1:nx+1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),6]   = np.reshape(inodes[1:nz+1,  1:ny+1, 1:nx+1],(nz)*(ny)*(nx),1) 
#    elements[0:(nx*ny*nz),7]   = np.reshape(inodes[1:nz+1,  1:ny+1, 0:-1],(nz)*(ny)*(nx),1)  
##    # edges
##    cnt = 0
##    # E0-------------------------------------------
##    E0  = inodes[0,:,:]
##    E00 = np.reshape(E0[0:ny,0:nx],nx*ny)
##    E01 = np.reshape(E0[0:ny,1:nx+1],nx*ny)
##    E02 = np.reshape(E0[1:ny+1,1:nx+1],nx*ny)
##    E03 = np.reshape(E0[1:ny+1,0:nx],nx*ny)
##    ind=np.arange(0,nx*ny)+cnt
##    cnt = cnt + nx*ny
##    edges[ind,0] = E00
##    edges[ind,1] = E01
##    edges[ind,2] = E02
##    edges[ind,3] = E03
##    edges[ind,4] = 0
##    # E1-------------------------------------------
##    E1  = inodes[-1,:,:] 
##    E10 = np.reshape(E1[0:ny,0:nx],nx*ny)
##    E11 = np.reshape(E1[0:ny,1:nx+1],nx*ny)
##    E12 = np.reshape(E1[1:ny+1,1:nx+1],nx*ny)
##    E13 = np.reshape(E1[1:ny+1,0:nx],nx*ny) 
##    ind=np.arange(0,(nx*ny))+cnt
##    cnt = cnt + nx*ny
##    edges[ind,0] = E10
##    edges[ind,1] = E13
##    edges[ind,2] = E12
##    edges[ind,3] = E11
##    edges[ind,4] = 1
##    # E2-------------------------------------------
##    E2  = inodes[:,0,:] 
##    E20 = np.reshape(E2[1:nz+1,0:nx],nx*nz)
##    E21 = np.reshape(E2[0:nz,0:nx],nx*nz)
##    E22 = np.reshape(E2[0:nz,1:nx+1],nx*nz)
##    E23 = np.reshape(E2[1:nz+1,1:nx+1],nx*nz) 
##    ind=np.arange(0,nx*nz)+cnt
##    cnt = cnt + nx*nz
##    edges[ind,0] = E20
##    edges[ind,1] = E23
##    edges[ind,2] = E22
##    edges[ind,3] = E21       
##    edges[ind,4] = 2     
##    # E3-------------------------------------------
##    E3  = inodes[:,-1,:] 
##    E30 = np.reshape(E3[0:nz,0:nx],nx*nz)
##    E31 = np.reshape(E3[1:nz+1,0:nx],nx*nz)
##    E32 = np.reshape(E3[1:nz+1,1:nx+1],nx*nz)
##    E33 = np.reshape(E3[0:nz,1:nx+1],nx*nz) 
##    ind=np.arange(0,nx*nz)+cnt
##    cnt = cnt + nx*nz
##    edges[ind,0] = E30
##    edges[ind,1] = E33
##    edges[ind,2] = E32
##    edges[ind,3] = E31
##    edges[ind,4] = 3
##    # E4-------------------------------------------
##    E4  = inodes[:,:,0] 
##    E40 = np.reshape(E4[0:nz,0:ny],ny*nz)
##    E41 = np.reshape(E4[0:nz,1:ny+1],ny*nz)
##    E42 = np.reshape(E4[1:nz+1,1:ny+1],ny*nz)
##    E43 = np.reshape(E4[1:nz+1,0:ny],ny*nz) 
##    ind=np.arange(0,ny*nz)+cnt
##    cnt = cnt + ny*nz
##    edges[ind,0] = E40
##    edges[ind,1] = E41
##    edges[ind,2] = E42
##    edges[ind,3] = E43         
##    edges[ind,4] = 4   
##    # E5-------------------------------------------
##    E5  = inodes[:,:,-1] 
##    E50 = np.reshape(E5[1:nz+1,0:ny],ny*nz)
##    E51 = np.reshape(E5[1:nz+1,1:ny+1],ny*nz)
##    E52 = np.reshape(E5[0:nz,1:ny+1],ny*nz)
##    E53 = np.reshape(E5[0:nz,0:ny],ny*nz) 
##    ind=np.arange(0,ny*nz)+cnt
##    cnt = cnt + ny*nz
##    edges[ind,0] = E50
##    edges[ind,1] = E51
##    edges[ind,2] = E52
##    edges[ind,3] = E53
##    edges[ind,4] = 5             
#
#    if Nx*Ny*Nz>1:
#        #a[ix_([1,3,4],[0,2])]
#        X=np.zeros((nx,ny,nz),np.int32)
#        c=0
#        for i in range(Nx):
#            for j in range(Ny):
#                for k in range(Nz):
#                    nNx = int(nx/Nx) 
#                    nNy = int(ny/Ny)
#                    nNz = int(nz/Nz)
#                    in1 = np.arange((i)*nNx,(i+1)*nNx)
#                    in2 = np.arange((j)*nNy,(j+1)*nNy)
#                    in3 = np.arange((k)*nNz,(k+1)*nNz)
#                    X[np.ix_(in1,in2,in3)]=c
#                    c=c+1
#
#        subs0 = X.reshape(nx*ny*nz,order='F')
#    else:
#        subs0 = np.zeros(1)
#
#    return {'coordinates':coordinates,'elements': elements, 'edges':edges,'decomp':subs0}
#
#    #################################################################
#    #                                                               #
#    #                         A z                                   #  
#    #                         |                                     #    
#    #                         |            L4                       #
#    #                         |_ _ _ _ _ _ _                        #
#    #                        /     L1      /|                       #  
#    #                       /_ _ _ _ _ _  / |                       #
#    #                      |      |      |  |                       #  
#    #                    L2|      |      |L3|                       #   
#    #                      |_ _ _ |_ _ _ |  |                       #  
#    #                      |    L5|      |  |-------> y             #
#    #                      |      |      | /                        #
#    #                      |_ _ _ |_ _ _ |/                         #
#    #                     /                                         #  
#    #                    /       L0                                 #
#    #                   /                                           #  
#    #                  v  x                                         #
#    #                                                               #
#    #################################################################   
#
#
#
