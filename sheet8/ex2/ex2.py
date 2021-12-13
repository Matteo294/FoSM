import numpy as np
from matplotlib import pyplot as plt 

N1 = 9
N2 = 5
N3 = 3
N4 = 2
L = 1

grid1 = np.linspace(0, L, N1)
grid2 = np.linspace(0, L, N2)
grid3 = np.linspace(0, L, N3)
grid4 = np.linspace(0, L, N4)

I_2h_h = np.array([1/2, 1, 1/2])
I_h_2h = np.array([1/4, 1/2, 1/4])

R1 = np.zeros((N2, N1)) # from grid 1 to grid 1
R2 = np.zeros((N3, N2)) # from grid 2 to grid 3
R3 = np.zeros((N4, N3)) # from grid3 to grid 4

for i in range(1, N2-1):
    R1[i,2*i-1] = 1
    R1[i,2*i] = 2
    R1[i,2*i+1] = 1
R1[0][0] = 2
R1[0][1] = 1
R1[N2-1][2*(N2-1)] = 2
R1[N2-1][2*(N2-1)-1] = 1
for i in range(1, N3-1):
    R2[i,2*i-1] = 1
    R2[i,2*i] = 2
    R2[i,2*i+1] = 1
R2[0][0] = 2
R2[0][1] = 1
R2[N3-1][2*(N3-1)] = 2
R2[N3-1][2*(N3-1)-1] = 1
for i in range(1, N4-1):
    R3[i,2*i-1] = 1
    R3[i,2*i] = 2
    R3[i,2*i+1] = 1 
R3[0][0] = 2
R3[0][1] = 1
R3[N4-1][2*(N4-1)] = 2
R3[N4-1][2*(N4-1)-1] = 1

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

A2 = np.einsum("ij,jk,kl->il", R1, A1, P1)
A3 = np.einsum("ij,jk,kl->il", R2, A2, P2)
A4 = np.einsum("ij,jk,kl->il", R3, A3, P3)

print(A1)
print(A2)
print(A3)
print(A4)