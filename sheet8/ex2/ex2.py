import numpy as np
from matplotlib import pyplot as plt 

N1 = 9
N2 = 5
N3 = 3
N4 = 2
L = 1
D = 1
h = L/N1
eps = 1
nrelax = 5
ncycles = 1

def Jacobi_Iteration(x, Dinv, L, U, b):
    #print(x.shape, Dinv.shape, L.shape, U.shape, b.shape)
    return Dinv*b + Dinv*np.einsum("ij,j", L+U, x)

# Matrix to move from 2h to h
def create_restr_mat(N):
    assert(N%2 == 1)
    Nrestr = int((N-1)/2) + 1
    R = np.zeros((Nrestr, N))
    for i in range(1, Nrestr-1):
        R[i,2*i-1] = 1
        R[i,2*i] = 2
        R[i,2*i+1] = 1
    R[0][0] = 2
    R[0][1] = 1
    R[Nrestr-1][2*(Nrestr-1)] = 2
    R[Nrestr-1][2*(Nrestr-1)-1] = 1
    return R

# Returns M1*A*M2
def restrict_prolong_matrix(A, M1, M2):
    return np.einsum("ij,jk,kl->il", M1, A, M2)
# Returns M*v
def restrict_prolong_vector(M, v):
    return np.einsum("ij,j", M, v)

# Returns D^-1, L, U given matrix A
def Dinv_L_U_decomposition(A):
    D_inv = np.zeros(len(A))
    L = np.zeros((len(A), len(A)))
    U = np.zeros((len(A), len(A)))
    for i in range(len(A)):
        for j in range(len(A)):
            if i < j:
                L[i,j] = -A[i,j]
            elif i > j:
                U[i,j] = -A[i,j]
            elif i == j:
                D_inv[i] = 1/A[i,i]
    return D_inv, L, U

grid1 = np.linspace(0, L, N1)
grid2 = np.linspace(0, L, N2)
grid3 = np.linspace(0, L, N3)
grid4 = np.linspace(0, L, N4)

R1 = create_restr_mat(N1)
R2 = create_restr_mat(N2)
R3 = create_restr_mat(N3)
P1 = 2*np.transpose(R1)
P2 = 2*np.transpose(R2)
P3 = 2*np.transpose(R3)

A1 = np.zeros((N1, N1))
for i in range(1, N1-1):
    A1[i,i-1] = 1
    A1[i,i] = -2
    A1[i,i+1] = 1
A1[0][0] = 1
A1[N1-1][N1-1] = 1

A2 = restrict_prolong_matrix(A1, R1, P1)
A3 = restrict_prolong_matrix(A2, R2, P2)
A4 = restrict_prolong_matrix(A3, R3, P3)

b1 = np.full(N1, -h**2/D/eps)
b1[0] = 1
b1[N1-1] = 1
b2 = restrict_prolong_vector(R1, b1)
b3 = restrict_prolong_vector(R2, b2)
b4 = restrict_prolong_vector(R3, b3)

D_inv1, L1, U1 = Dinv_L_U_decomposition(A1)
D_inv2, L2, U2 = Dinv_L_U_decomposition(A2)
D_inv3, L3, U3 = Dinv_L_U_decomposition(A3)
D_inv4, L4, U4 = Dinv_L_U_decomposition(A4)

x1 = np.zeros(N1)

for n in range(ncycles):
    print("Start v-cycle", n+1, ", Nrelax:", nrelax, ", Ncycles:", ncycles)
    for _ in range(nrelax):
        x1 = Jacobi_Iteration(x1, D_inv1, L1, U1, b1)
    x2 = restrict_prolong_vector(R1, x1)
    for _ in range(nrelax):
        x2 = Jacobi_Iteration(x2, D_inv2, L2, U2, b2)
    x3 = restrict_prolong_vector(R2, x2)
    for _ in range(nrelax):
        x3 = Jacobi_Iteration(x3, D_inv3, L3, U3, b3)
    x4 = restrict_prolong_vector(R3, x3)
    for _ in range(nrelax):
        x4 = Jacobi_Iteration(x4, D_inv4, L4, U4, b4)
    x3 = restrict_prolong_vector(P3, x4)
    for _ in range(nrelax):
        x3 = Jacobi_Iteration(x3, D_inv3, L3, U3, b3)
    x2 = restrict_prolong_vector(P2, x3)
    for _ in range(nrelax):
        x2 = Jacobi_Iteration(x2, D_inv2, L2, U2, b2)
    x1 = restrict_prolong_vector(P1, x2)
    for _ in range(nrelax):
        x1 = Jacobi_Iteration(x1, D_inv1, L1, U1, b1)
    print("Mean squared error", sum(np.einsum("ij,j", A1, x1) - b1)**2)
