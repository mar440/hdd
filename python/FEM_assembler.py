from __future__ import division
from __future__ import print_function

import numpy as np

def stima3e8(inputs):

    coordinates = inputs['coordinatesLoc']
    e           = inputs['Ex']
    mu          = inputs['mu']
    volF        = inputs['volF']

    
#    u           = inputs['u']
#    du          = inputs['du']
#    material    = inputs['material']
#    stressIn    = inputs['stressIn']
#    stressDin   = inputs['stressDin']
#    epelIn      = inputs['epelIn']
#    epplIn      = inputs['epplIn']
#    epeqIn      = inputs['epeqIn']
#    plworkIn    = inputs['plworkIn']
#    statevIn    = inputs['statevIn']
#    YvIn        = inputs['YvIn']
    
    nP  = coordinates.shape[0]
    dim = 3
    neq = nP * dim 


    d = np.zeros((6,6))
    const   = e/((1.+mu)*(1.-2.*mu));
    mu2     = (1-mu)
    mu3     = (0.5-mu)
    d[0,]=np.array([mu2,    mu,     mu,     0.,     0.,     0.])*const
    d[1,]=np.array([mu,     mu2,    mu,     0.,     0.,     0.])*const
    d[2,]=np.array([mu,     mu,     mu2,    0.,     0.,     0.])*const
    d[3,]=np.array([0.,     0.,     0.,     mu3,    0.,     0.])*const
    d[4,]=np.array([0.,     0.,     0.,     0.,     mu3,    0.])*const
    d[5,]=np.array([0.,     0.,     0.,     0.,     0.,     mu3])*const

    pt = (1./np.sqrt(3))
    R = np.array([[-pt,-pt,-pt],\
            [-pt,-pt, pt],\
            [-pt, pt,-pt],\
            [-pt, pt, pt],\
            [ pt,-pt,-pt],\
            [ pt,-pt, pt],\
            [ pt, pt,-pt],\
            [ pt, pt, pt]])

    B = np.zeros((6,neq))
    NN = np.zeros((neq,3))
    Ke = np.zeros((neq,neq))
    fe = np.zeros(Ke.shape[0])

    indG = np.arange(0,6)
    indStatev = np.arange(0,6*np.shape(d)[0])


    for i in range(8):
        r = R[i,0]; s = R[i,1]; t = R[i,2]

        dNr=0.125*np.array([-(1.-s)*(1.-t), (1.-s)*(1.-t), (1.+s)*(1.-t),-(1.+s)*(1.-t),\
                -(1.-s)*(1.+t), (1.-s)*(1.+t), (1.+s)*(1.+t),-(1.+s)*(1.+t)])
        dNs=0.125*np.array([-(1.-r)*(1.-t),-(1.+r)*(1.-t), (1.+r)*(1.-t), (1.-r)*(1.-t),\
                -(1.-r)*(1.+t),-(1.+r)*(1.+t), (1.+r)*(1.+t), (1.-r)*(1.+t)]);
        dNt=0.125*np.array([-(1.-r)*(1.-s),-(1.+r)*(1.-s),-(1.+r)*(1.+s),-(1.-r)*(1.+s),\
                (1.-r)*(1.-s), (1.+r)*(1.-s), (1.+r)*(1.+s), (1.-r)*(1.+s)]);

        N = np.array([(1.-r)*(1.-s)*(1.-t),(r+1.)*(1.-s)*(1.-t),(r+1.)*(s+1.)*(1.-t),(1.-r)*(s+1.)*(1.-t),\
                (1.-r)*(1.-s)*(t+1.),(r+1.)*(1.-s)*(t+1.),(r+1.)*(s+1.)*(t+1.),(1.-r)*(s+1.)*(t+1.)]);                     

        J = np.dot(np.array([dNr,dNs,dNt]),coordinates)
        iJ = np.linalg.inv(J)
        dJ = np.abs(np.linalg.det(J))

        dNx = iJ[0,0]*dNr + iJ[0,1]*dNs + iJ[0,2]*dNt
        dNy = iJ[1,0]*dNr + iJ[1,1]*dNs + iJ[1,2]*dNt
        dNz = iJ[2,0]*dNr + iJ[2,1]*dNs + iJ[2,2]*dNt 

        ZerosMat = np.zeros(8)
        B[0,0:8]   = dNx
        B[1,8:16]  = dNy
        B[2,16:24] = dNz
        B[3,0:8]   = dNy;  B[3,8:16]  = dNx
        B[4,8:16]  = dNz;  B[4,16:24] = dNy
        B[5,0:8]   = dNz;  B[5,16:24] = dNx


        Ke += np.dot(B.transpose(),np.dot(d,B))*dJ
        fe[0*nP:1*nP] += N * volF[0]*dJ 
        fe[1*nP:2*nP] += N * volF[1]*dJ 
        fe[2*nP:3*nP] += N * volF[2]*dJ 

        indG        = indG + 6
        indStatev   = indStatev + 6*d.shape[0]
        flagElast   = True

    return {'Ke':Ke,'fe': fe}





