import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

class Density_map:
    def __init__(self, H, K, m=1):
        self.H = H # Grid goes from [-H,H] in both dimensions
        self.K = K # number of grid cells per dimension
        self.h = 2*H/K # cell size per dimension
        self.m = m # mass per particle, assume all particles same mass
        self.particles = [] # list of particle coordinates
        self.W0 = np.zeros((K,K)) # zeroth order weight matrix, each entry corresponds to a cell
        self.W1 = np.copy(self.W0) # same for 1st order
        self.W2 = np.copy(self.W0) # 2nd order

        self.rho0 = np.copy(self.W0) # also initiate rho matrices so that they can be later called by grid.rho
        self.rho1 = np.copy(self.W0)
        self.rho2 = np.copy(self.W0)

    class Particle:
        def __init__(self,X,grid):
            self.grid = grid # needed to be able to call Particle.grid.get_indices_k_l
            self.x, self.y = X[0], X[1] # cartesian coordinates
            self.X = X # vector with x and y cartesian coords
            self.k, self.l = self.grid.get_indices_k_l(self.X)[0] # theoretical k l indices
            self.ly = self.grid.get_indices_k_l(self.X)[1] # y index ly because in matrix, we go from up to down but in physical picture from down to up

    def get_indices_k_l(self,X):
        k = int(np.floor((X[0]+self.H)/self.h)) 
        l_index = -int(np.floor((X[1]+self.H)/self.h)) - 1 # y index l for indexing bc matrix goes other direction than physical picture
        l_theory = int(np.floor((X[1]+self.H)/self.h)) # l index as in theory starting to count from down to up
        #print(k,l_theory)
        return (k, l_theory), l_index

    def top_hat(self, x): # top hat function // 1 for |x| <= 1/2, 0 else
        if np.abs(x) <= 1/2:
            return 1
        else:
            return 0

    def get_eps(self,particle): # calculate epsx and y: (X_i - x_{k-1/2})/h
        ex = (particle.x - self.H*(2*particle.k/self.K - 1))/self.h
        ey = (particle.y - self.H*(2*particle.l/self.K - 1))/self.h
        return ex, ey

    def add_particle(self, X): # add particle to grid.particles list, has properties cartesian coordinates and k l and ly indices // updates weight matrix and thus rho for each newly added particle
        self.particles.append(self.Particle(X,self))
        k, l = self.particles[-1].k, self.particles[-1].ly
    
#        if k == 0 or k == K-1 or l == 0 or l == K-1: # for simplicity just discard particles laying outside in plot_map.py
#           maybe add error message and exit this function?   

        self.W0[l,k] += 1 # 0th order, simply update k l cell where particle sits
        
        ex, ey = self.get_eps(self.particles[-1]) 
        
        # 1st order, touches 4 cells in total, which 4 out of 9 possible is determined by top hat of ex and ey
        self.W1[l,k] += ex*ey
        self.W1[l+1,k] += ex*(1-ey)*(1-self.top_hat(ey))
        self.W1[l-1,k] += ex*(1-ey)*self.top_hat(ey)
        self.W1[l,k+1] += (1-ex)*(1-self.top_hat(ex))*ey
        self.W1[l,k-1] += (1-ex)*self.top_hat(ex)*ey
        self.W1[l+1,k+1] += (1-ex)*(1-self.top_hat(ex))*(1-ey)*(1-self.top_hat(ey))
        self.W1[l+1,k-1] += (1-ex)*self.top_hat(ex)*(1-ey)*(1-self.top_hat(ey))
        self.W1[l-1,k+1] += (1-ex)*(1-self.top_hat(ex))*(1-ey)*self.top_hat(ey)
        self.W1[l-1,k-1] += (1-ex)*self.top_hat(ex)*(1-ey)*self.top_hat(ey)

        # 2nd order, touches all 9 cells (8 neighbours), see geometrical derivation of the expressions
        self.W2[l,k] += (0.5 + ex - ex**2)*(0.5 + ey - ey**2)
        self.W2[l+1,k] += (0.5 + ex - ex**2)*(0.5*ey**2)
        self.W2[l-1,k] += (0.5 + ex - ex**2)*((1-ey)**2/2)
        self.W2[l,k+1] += (0.5*ex**2)*(0.5 + ey - ey**2)
        self.W2[l,k-1] += ((1-ex)**2/2)*(0.5 + ey - ey**2)
        self.W2[l+1,k+1] += (0.5*ex**2)*(0.5*ey**2) 
        self.W2[l+1,k-1] += ((1-ex)**2/2)*(0.5*ey**2)
        self.W2[l-1,k+1] += (0.5*ex**2)*((1-ey)**2/2)
        self.W2[l-1,k-1] += ((1-ex)**2/2)*((1-ey)**2/2)               

        self.rho0 = self.m*self.W0/self.h**2
        self.rho1 = self.m*self.W1/self.h**2
        self.rho2 = self.m*self.W2/self.h**2
        
    def show_density_map(self, rho): # shows rho as a heatmap, colorbar is the density per cell
        ax = sns.heatmap(rho,vmin=0.0,cmap='inferno')
        ticks = np.arange(0,self.K+1,1) # set new x and y ticks corresponding to the physical dimensions going from -H to H
        names = np.arange(-self.H,self.H+1,1)
        ax.set_xticks(ticks=ticks)
        ax.set_yticks(ticks=ticks)
        ax.set_yticklabels(labels=-names) # - bc should start from below not from above
        ax.set_xticklabels(labels=names)
        plt.grid()
        plt.show()
