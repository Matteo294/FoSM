import numpy as np
from Density_map import *

# 2)

H = 10.0
K = 20

def ex2():

    grid = Density_map(H,K)
    x, y = np.random.uniform(low=-grid.H + grid.h, high = grid.H - grid.h, size=2)
    X = np.array((x,y))

    print(X)

    grid.add_particle(X)

    grid.show_density_map(grid.rho0)
    grid.show_density_map(grid.rho1)
    grid.show_density_map(grid.rho2)

def ex34(N):
    
    grid = Density_map(H,K)
    for n in range(N):
        x, y = np.random.normal(0.0,3.0,2)
        if np.abs(x) < grid.H - grid.h and np.abs(y) < grid.H - grid.h:
            X = np.array((x,y))
            grid.add_particle(X)

    grid.show_density_map(grid.rho0)
    # integral over rho
    M0 = grid.rho0.sum()
    print(f'0th order: M = {M0}')
    grid.show_density_map(grid.rho1)
    M1 = grid.rho1.sum()
    print(f'1st order: M = {M1}')
    grid.show_density_map(grid.rho2) 
    M2 = grid.rho2.sum()
    print(f'2nd order: M = {M2}')



if __name__ == "__main__":
    ex2()
    ex34(1000)
