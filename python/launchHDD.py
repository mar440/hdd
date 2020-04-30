import numpy as np
import pylab as plt
import readCOO
from scipy.sparse.linalg import spsolve
import singularCholesky as cholesky
import json


with open("../hddConf.json", "r") as read_file: 
    hddConf = json.load(read_file) 

meshOpt = hddConf['builtInTestCaseSetting']['mesh']
solverOpt = hddConf['solver']


nSx = meshOpt['numberOfSubdomains_x']
nSy = meshOpt['numberOfSubdomains_y']

precond = solverOpt['preconditioner']

nSUB = nSx * nSy
offset = 1 # COO sparse matrix market format begins with '1'
makeSymmetric=False


nDOFs = []


def readTxtFiles(pathBegin):
    u=[]
    K=[]
    f=[]
    B=[]
    R=[]
    BCOO=[]
    dataB=[]
    nullPivots=[]

    for iSub in range(nSUB):

        path = pathBegin+"/solution_"+str(iSub)+".txt"
        _u = np.loadtxt(path).astype(float)
        u.append(_u)

        path = pathBegin+"/K_"+str(iSub)+".txt"
        _K=readCOO.load_matrix_basic(path,offset,makeSymmetric)
        K.append(_K)

        path = pathBegin+"/R_"+str(iSub)+".txt"
        _R = np.loadtxt(path)
        R.append(_R)

        path = pathBegin+"/B_"+str(iSub)+".txt"
        _B = np.loadtxt(path).astype(float)
        dataB.append(_B)

        path = pathBegin+"/f_"+str(iSub)+".txt"
        _f = np.loadtxt(path)
        f.append(_f)

        path = pathBegin+"/nullPivots_"+str(iSub)+".txt"
        _nlpv = np.loadtxt(path)
        nullPivots.append(_nlpv)

        nDOFs.append(_K.shape[0])
        BCOO.append(np.zeros((0,3),dtype=int))
        B.append(np.zeros((0,0)))

    return u, K, R, dataB, f, nullPivots, B, BCOO




def createBooleanMatrix():

    cntLambda = 0
    scaling = np.zeros(0)
    for iSub in range(nSUB):
        dataB_i =  dataB[iSub]
        neighbours = np.unique(dataB_i[:,0]).astype(int)
        for cnt, jSub in enumerate(neighbours):
            iRows_subB_i = dataB_i[:,0] == jSub 
            subSetB_i = dataB_i[iRows_subB_i,:]
            if (subSetB_i[0,-1] > 0):
                data_B_j = dataB[jSub]
                iRows_subB_j = data_B_j[:,0] == iSub 
                subSetB_j = data_B_j[iRows_subB_j,:]
                for rw in range(subSetB_i.shape[0]):
                    BCOO[iSub] = np.vstack((BCOO[iSub], np.array([cntLambda, subSetB_i[rw,1], subSetB_i[rw,2]])))
                    BCOO[jSub] = np.vstack((BCOO[jSub], np.array([cntLambda, subSetB_j[rw,1], subSetB_j[rw,2]])))
                    scaling = np.concatenate((scaling,[1.0 / subSetB_i[rw,-1]]))
                    cntLambda += 1

    for iSub in range(nSUB):
        B[iSub] = np.zeros((cntLambda,nDOFs[iSub]))
        for rwB in range(BCOO[iSub].shape[0]):
            I = BCOO[iSub][rwB,0].astype(int)
            J = BCOO[iSub][rwB,1].astype(int)
            V = BCOO[iSub][rwB,2]
            B[iSub][I,J] = V
    return B, scaling, cntLambda



def createDirichletPrecond(listOfK):
    DirichletPrecond = []
    for iSub in range(nSUB):
        inDOFs = listOfK[iSub].shape[0]
        indGamma = np.unique(BCOO[iSub][:,1].astype(int))
        indIner = np.arange(inDOFs)
        indIner = np.delete(indIner,indGamma)
        KBB = listOfK[iSub][indGamma,:][:,indGamma]     # [np.ix_(indGamma,indGamma)]
        KII = listOfK[iSub][indIner ,:][:,indIner]      # [np.ix_(indIner,indIner)]
        KIB = listOfK[iSub][indIner ,:][:,indGamma]     # [np.ix_(indIner,indGamma)]
        KBI = listOfK[iSub][indGamma,:][:,indIner]      # [np.ix_(indGamma,indIner)]
        SBB = KBB - KBI.dot(spsolve(KII,KIB))
        DirichletPrecond.append(SBB)
    return DirichletPrecond


def nonSingularStiffMats():

    Kreg  = []
    for iSub in range(nSUB):
        Kreg.append(K[iSub].copy())

    for iSub in range(nSUB):
        nlPv = nullPivots[iSub]
        if (nlPv.shape[0] > 0):
            for pivot in nlPv:
                Kreg[iSub][pivot,pivot] *= 2
    return Kreg

def createDualObjects():

    F = np.zeros((cntLambda,cntLambda))
    d = np.zeros((cntLambda))
    G = np.zeros((cntLambda,0))
    e = np.zeros((0))

    for iSub in range(nSUB):
        F += B[iSub].dot(spsolve(Kreg[iSub],B[iSub].T))
        Kplusf = spsolve(Kreg[iSub],f[iSub])

        d += B[iSub].dot(Kplusf)
        if (R[iSub].shape[0] > 0):
            Gi = -B[iSub].dot(R[iSub])
            G = np.hstack((G,Gi.copy()))
            ei = -R[iSub].T.dot(f[iSub])
            e = np.concatenate((e,ei)) 

    return F,d,G,e

############################################################
############################################################
############################################################


pathBegin = "/home/mar440/WorkSpace/hfls/mfeti/build/"
#pathBegin = "../build/"
u, K, R, dataB, f, nullPivots, B, BCOO = readTxtFiles(pathBegin)
B,scaling, cntLambda  = createBooleanMatrix()
DirichletPrecond = createDirichletPrecond(K)
Kreg = nonSingularStiffMats()
F,d,G,e = createDualObjects()

invGtG = np.linalg.inv(G.T.dot(G))

def Proj(x):
    alpha = -invGtG.dot(G.T.dot(x))
    Px = x + G.dot(alpha)
    return Px, alpha

def Precond(lam_inp, whichType, indSubdomain=-1):

    B_scaling = np.diag(scaling)


    if indSubdomain == -1:
        indexXlist = np.arange(nSUB).tolist()
    else:
        indexXlist = [indSubdomain]

    Zi = np.zeros(lam_inp.shape)
    for iSub_ in indexXlist:
        Bt_lam = B[iSub_].T.dot(B_scaling.dot(lam_inp))

        if (whichType == 'Dirichlet'):
            Sx = np.zeros(B[iSub_].shape[1])
            indGamma = np.unique(BCOO[iSub_][:,1].astype(int))
            Sx[indGamma] = DirichletPrecond[iSub_].dot(Bt_lam[indGamma])

            SBtX1 = K[iSub_].dot(Bt_lam)
#            for ii in indGamma:
#                print("aaa--- ", SBtX1[ii], " ", Sx[ii])


        elif (whichType == 'Lumped'):
            Sx = K[iSub_].dot(Bt_lam)
        else:
            print("errrrror in precond")
            return []

        Zi += B_scaling.dot(B[iSub_].dot(Sx))
#        Zi1 = B_scaling.dot(B[iSub_].dot(K[iSub_].dot(Bt_lam)))
#        for ii in range(Zi.shape[0]):
#            print("aaa--- ", Zi[ii], " ", Zi1[ii])

    return Zi





def singsolve(_A,_b):

    nA = _A.shape[0]

    dd = np.diag(np.sqrt(1.0 / np.diag(_A)))
    scaledA = dd.dot(_A.dot(dd))
    scaledRHS = dd.dot(_b)

    if True:

        L, P, rank = cholesky.ldl_factorP(scaledA)
#        print("rank\n",rank)
        if False:#rank < scaledA.shape[0]:
            print(dd)
            print("scaledA:\n",scaledA)
            print("A:\n",_A)
            print("rank: ", rank," nA: ", scaledA.shape[0]);

#        U,s,V = np.linalg.svd(scaledA)
#        for ss in s: print(ss.round(5),end=" ")
#        print(" ")

        Lr = L[:rank,:rank]
        Prhs = P.dot(scaledRHS)
        if (len(Prhs.shape) == 1):
            Prhs = Prhs[:rank]
        else:
            Prhs = Prhs[:rank,:]


#        print("rank: ", rank)
        y = np.linalg.solve(Lr,Prhs)
        out = P.T.dot(np.linalg.solve(Lr.T,y))


#        permutation =  p[:rank]
#
#        if len(scaledRHS.shape) == 1:
#            br = scaledRHS[permutation]
#        else:
#            br = scaledRHS[permutation,:]
#
#        Lreduced = L[:rank,:rank]
#        out_p = np.linalg.solve(Lreduced.T,np.linalg.solve(Lreduced,br))





        return  dd.dot(out)
    else:



        U,s,V = np.linalg.svd(scaledA)

        indx = s<1e-6
        s[indx] = 1.0
        invs = 1.0 / s
        invs[indx] = 0
        iS = np.diag(invs)

        return dd.dot(V.T.dot(iS.dot(U.T.dot(scaledRHS))))



#########################################################################################################
def PCPG(solverSetting):

    PrecondType = solverSetting[0]
    niter = solverSetting[1]
    tol = solverSetting[2]

#    S = np.zeros((cntLambda, cntLambda))
    B_scaling = np.diag(scaling)
#    for iSub in range(nSUB):
#        S += B[iSub].dot(K[iSub].dot(B[iSub].T))
        #S += B[iSub].dot(B[iSub].T)
#    S = B_scaling.dot(S.dot(B_scaling))

    lambda0 = G.dot(invGtG.dot(e))
    g = F.dot(lambda0) - d
    Pg, alpha = Proj(g)
    normPg0 = np.linalg.norm(Pg)
    Z = Precond(Pg,PrecondType)
#    Z = S.dot(Pg)
    W,beta = Proj(Z)
    Lambda = lambda0.copy()
    gamma = Z.T.dot(Pg) 
    sqrt_gamma0 = np.sqrt(gamma.copy())
    vectorNorm = []

    print("|| gamma0 || = " ,sqrt_gamma0, ", tol = " , tol)


    for i in range(niter):

        vectorNorm.append(np.sqrt(gamma) /sqrt_gamma0)
        print("i: ", i+1, " ||ztPg|| = " , np.sqrt(gamma) /sqrt_gamma0, end=" ")
        print(" ||gtPg|| = " , np.linalg.norm(Pg) /normPg0)
        if ( (np.sqrt(gamma) /sqrt_gamma0) < tol):
            break

        Q = F.dot(W)
        Delta = Q.T.dot(W)
        rho = -gamma / Delta
        Lambda += W * rho
        g += Q.dot(rho)
        Pg, alpha = Proj(g)
        SPg = Precond(Pg,PrecondType)
#        SPg = S.dot(Pg)
        Z, _beta = Proj(SPg)

        gamma_prev = gamma

        gamma = Z.T.dot(Pg) 

        Pg_prv = Pg.copy()
        W_prv = W.copy()

        SPg = Precond(Pg,PrecondType)
        PZ, _beta = Proj(SPg)
        W = PZ + W_prv.dot(gamma / gamma_prev)

    upy = []
    offsetAlpha = 0
    for iSub in range(nSUB):
        defect = 0
        if R[iSub].shape[0]:
            defect = R[iSub].shape[1]
        alpha_i = alpha[offsetAlpha:(offsetAlpha + defect)]
        offsetAlpha += defect
        _u = spsolve(Kreg[iSub],f[iSub] - B[iSub].T.dot(Lambda)) + R[iSub].dot(alpha_i)
        upy.append(_u)

    print("F*l - d + Gt*a = ",np.linalg.norm(F.dot(Lambda) - d + G.dot(alpha)))

    return Lambda, alpha, upy, np.asarray(vectorNorm)
#
#########################################################################################################
def MPCPG(solverSetting):

    PrecondType = solverSetting[0]
    niter = solverSetting[1]
    tol = solverSetting[2]

    lambda0 = G.dot(invGtG.dot(e))
    g = F.dot(lambda0) - d

    Z = np.zeros((cntLambda,nSUB))
#    B_scaling = np.diag(scaling)

    Pg, alpha = Proj(g)
    normPg0 = np.linalg.norm(Pg)
    for iSub in range(nSUB):
        Zi = Precond(Pg,PrecondType,iSub)
        Z[:,iSub] = Zi

    W, _beta = Proj(Z)
    Lambda = lambda0.copy()

    listQ = []
    listW = [W]

    normGammaPg0 = 1

    vectorNorm = []

    for Iter in range(niter):

        W_i = listW[-1]
        Q = F.dot(W_i)
        Delta_i = Q.T.dot(W_i)

        L1,N1,rank1 = cholesky.ldl_factorP(Delta_i)
        iL1 = np.zeros((rank1,L1.shape[0]))
        iL1[:rank1,:rank1] = np.linalg.inv(L1[:rank1,:rank1])
        if (rank1 < Delta_i.shape[0]):
            print("defect: ",Delta_i.shape[0] - rank1)

        W_i = W_i.dot(N1.T.dot(iL1.T))
        listW[-1] = W_i

        Q = Q.dot(N1.T.dot(iL1.T))
        listQ.append(Q)

        gamma = W_i.T.dot(Pg)
        rho = (-1) * gamma

        if Iter == 0:
            normGammaPg0 = np.sum(gamma)
            normCG = 1;
        else:
            normCG = abs(np.sum(gamma))/normGammaPg0
        print("i: ", Iter + 1, " ||ZtPg|| = " , normCG )


        vectorNorm.append(normCG)

        if (normCG < tol):
            break


        Lambda += W_i.dot(rho)
        g += Q.dot(rho)
        Pg, alpha = Proj(g)
#        for iSub in range(nSUB):
#            Zi = B_scaling.dot(B[iSub].dot(K[iSub].dot(B[iSub].T)).dot(B_scaling.dot(Pg)))
#            Z[:,iSub] =Zi

        for iSub in range(nSUB):
            Zi = Precond(Pg,PrecondType,iSub)
            Z[:,iSub] = Zi



        PZ_i = Proj(Z)[0]
        W_next = PZ_i.copy()

        for inerIter in range(Iter+1):
            Wk = listW[-1 - inerIter]
            Qk = listQ[-1 - inerIter]
            Phi = -Qk.T.dot(W_next)
            W_next += Wk.dot(Phi)

        listW.append(W_next)

    upy = []
    offsetAlpha = 0
    for iSub in range(nSUB):
        defect = 0
        if R[iSub].shape[0]:
            defect = R[iSub].shape[1]
        alpha_i = alpha[offsetAlpha:(offsetAlpha + defect)]
        offsetAlpha += defect
        _u = spsolve(Kreg[iSub],f[iSub] - B[iSub].T.dot(Lambda)) + R[iSub].dot(alpha_i)
        upy.append(_u)

    print("F*l - d + Gt*a = ",np.linalg.norm(F.dot(Lambda) - d + G.dot(alpha)))
    return Lambda, alpha, upy, np.asarray(vectorNorm)

def normOfJumps(displacements, label0):
    sumBu = 0
    for iSub in range(nSUB):
        sumBu += B[iSub].dot(displacements[iSub])

    print("|| B * u || = ", np.linalg.norm(sumBu), "   ", label0)
#

PrecondType = "DIRICHLET"
niter = 50
tol = 1e-6

solverSetting = [solverOpt['preconditioner'],\
                solverOpt['maxNumbIter'],\
                solverOpt['stopingCriteria']]


lam1, alp1,upy1, normG1 = PCPG(solverSetting)
lam2, alp2,upy2, normG2 = MPCPG(solverSetting)
##
ix1 = np.arange(normG1.shape[0])
ix2 = np.arange(normG2.shape[0])
plt.semilogy(ix1,normG1,'.-',ix2,normG2,'.-')
plt.legend(['PCPG','MPCPG'])
plt.show()
#
normOfJumps(upy1, "upy1")
normOfJumps(upy2, "upy2")
normOfJumps(u, "uc++")

