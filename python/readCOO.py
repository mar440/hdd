import numpy as np
from scipy.sparse import csc_matrix

def load_matrix_basic(pathToFile,offset,makeSymmetric=True):

    tmp = np.loadtxt(pathToFile, skiprows=1)

    if (tmp.shape[0]==1):
        tmp = []
    else:
        n = np.int32(tmp[0,0])   
        m = np.int32(tmp[0,1])
        I = tmp[1::,0]-offset;
        J = tmp[1::,1]-offset;
        V = tmp[1::,2]
        if (makeSymmetric):
            logInd = J != I;
            Inew = np.concatenate((I,J[logInd]))
            J = np.concatenate((J,I[logInd]))
            I = Inew 
            V = np.concatenate((V,V[logInd]))
        tmp = csc_matrix((V,(I,J)),shape=(n,m)).tocsc()
    return tmp

