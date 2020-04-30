import numpy as np


def ldl_factorP(A):

    ##  Pt * L * Lt * P = A
    ##  L * Lt = P * A * Pt
    counter = 0
    dD = np.diag(np.sqrt(np.diag(A)))
    iD = np.diag(1.0 / np.sqrt(np.diag(A)))
    L,p, counter = __ldl_factorP__(iD.dot(A.dot(iD)),counter)

    P = np.eye(A.shape[0])
    P = P[p,:]
    dD = P.dot(dD.dot(P.T))
    L = dD.dot(L)
    return L, P, counter


def __ldl_factorP__(A,counter):

    nullPivotTolerance = 1e-14
    n = A.shape[0]


    perm = np.arange(n)
    L = np.zeros(A.shape)
    dA = np.diag(A)
    p = np.argmax(dA)
    perm[-1] = p
    perm[p] = n-1
    Ac = A[np.ix_(perm,perm)]

    if Ac[0,0] < 0 or np.abs(Ac[0,0]) < nullPivotTolerance:
        return L, perm, counter


    #print("Ac = ", Ac[0,0])
    l11 = np.sqrt(Ac[0,0])
    L[0,0] = l11

    if n > 1:
        l1 = Ac[1:,0] / l11
        A_ = Ac[1:,1:] - np.outer(l1,l1) 
        L[1:,1:], pout, counter = __ldl_factorP__(A_,counter)
        perm[1:] = perm[1:][pout] 
        L[1:,0] = l1[pout]

    counter += 1

    return L, perm, counter

if __name__ == '__main__':

    A = np.array([
        [ 4.06462389e+02,  9.65043990e+02, -1.51768915e+01, 1.41151993e+02],
        [ 9.65043990e+02,  1.99313978e+04,  1.41151993e+02, 1.16119375e+04],
        [-1.51768915e+01,  1.41151993e+02,  4.06462389e+02, 9.65043990e+02],
        [ 1.41151993e+02,  1.16119375e+04,  9.65043990e+02, 1.99313978e+04]])

#    n = 5
#    A =np.random.rand(n,n)
#    U, s, V = np.linalg.svd(A)
#    s[-3:]=0
#    A = U.dot(np.diag(s).dot(U.T))
    L, P, rank = ldl_factorP(A)
    LLt_A = L.dot(L.T) - P.dot(A.dot(P.T))
    del0 =np.linalg.norm(LLt_A)
    print("L*Lt - A(p,p) = ", del0)
    print("rank = ", rank)
