import numpy as np


#def ldl_factorP(A):
#
#    L,p,defect = __ldl_factorP__(A)
#    return L, p, defect


def ldl_factorP(A):

    n = A.shape[0]
    defect = 0

    perm = np.arange(n)
    L = np.zeros(A.shape)
    dA = np.diag(A)
    p = np.argmax(dA)
    perm[-1] = p
    perm[p] = n-1
    Ac = A[np.ix_(perm,perm)]

    if np.abs(Ac[0,0]) < 1e-10:
        return L, perm, Ac.shape[0]


    l11 = np.sqrt(Ac[0,0])
    L[0,0] = l11

    if n > 1:
        l1 = Ac[1:,0] / l11
        A_ = Ac[1:,1:] - np.outer(l1,l1) 
        L[1:,1:], pout, defect = ldl_factorP(A_)
        perm[1:] = perm[1:][pout] 
        L[1:,0] = l1[pout]


    return L, perm, defect

if __name__ == '__main__':

    A = np.array([
        [ 1.        ,  0.21216552,  0.17292542, -0.10657882],
        [ 0.21216552,  1.        , -0.10657882,  0.66966913],
        [ 0.17292542, -0.10657882,  1.        ,  0.21216552],
        [-0.10657882,  0.66966913,  0.21216552,  1.        ]])

    n = 30 
    A =np.random.rand(n,n)
    U, s, V = np.linalg.svd(A)
    s[-6:]=0
    A = U.dot(np.diag(s).dot(U.T))
    L, p, defect  = ldl_factorP(A)
    LLt_A = L.dot(L.T) - A[np.ix_(p,p)]
    del0 =np.linalg.norm(LLt_A)
    print("L*Lt - A(p,p) = ", del0)
    print("defect = ", defect)
