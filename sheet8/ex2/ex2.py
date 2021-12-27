import numpy as np

# multigrid method applied to heat diffusion 1D problem
def Jacobi_step(A,b,x):
    
    n = A.shape()[0]
    D = np.diag(A)*np.eye(n)
    LU = A - D
    
    return np.linalg.inv(D)@b - np.linalg.inv(D)@LU@x


def restrict_mat(n,nn):
    # restrict matrix A to next lower level n --> nn
    R = np.zeros((n,nn))
    R[0,0] = R[-1,-1] = 1
    
    # return R
    pass

def prolong_mat(A):

    # prolong matrix A to next higher level

    # P = 
    # return P
    pass


def Vcycle_Jacobi(A,b,x):

    # use V shaped method  // recursively call restrict_mat / prolong_mat
    # ...
    # use Jacobi iteration
    x_1 = Jacobi_step(A,b,x)
    
    res_1 = b - A@x_1

    n = A.shape()[0]
    
    

    # x =
    # return x
    pass





