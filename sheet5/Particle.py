from matplotlib import pyplot as plt 
import numpy as np

###################### Grid class ########################
# Creates a 2D grid from -H to H divided into N^2 cells
# The side length of one square is h
# add particles to the grid with the subclass Particle
##########################################################
class Grid:
    def __init__(self, H, N, m=1):
        self.H = H
        self.N = N
        self.h = 2*H/N
        self.m = m
        self.particles = []
        self.W0 = np.zeros((N, N)) # weights matrix zeroth order
        self.W1 = np.zeros((N, N)) # weights matrix first order
        self.W2 = np.zeros((N, N)) # weights matrix second order

    class Particle:
        def __init__(self, r, grid_instance):
            self.grid = grid_instance
            self.x = r[0] # x coord
            self.y = r[1] # y coord
            self.r = r # position vector
            self.k, self.l = self.grid.find_kl(self.r) # get indices of the closest cell
            self.xh = self.x/self.grid.h # x coordinate in units of h
            self.yh = self.y/self.grid.h # y coordinate in units of h
            self.rg = [self.xh, self.yh] # r in units of h

    # Adds a particle to the grid at position r
    def add_particle(self, r):
        self.particles.append(self.Particle(r, self))
        k = self.particles[-1].k
        l = self.particles[-1].l
        # Update weights for zeroth order
        self.W0[k,l] += 1
        # Update weights for the first order
        epsx, epsy = self.get_eps(self.particles[-1])
        self.W1[k,l] += epsx*epsy
        self.W1[k,l+1] += 1
        self.W1[k,l-1] += 1
        self.W1[k+1,l] += 1
        self.W1[k+1,l+1] += 1
        self.W1[k+1,l-1] += 1
        self.W1[k-1,l] += 1
        self.W1[k-1,l+1] += 1
        self.W1[k-1,l-1] += 1
        # Update weights for the second order

    # Get indices of the cell closest to position r
    def find_kl(self, r):
        k = int((self.H+r[0])/self.h)
        l = int((self.H+r[1])/self.h)
        return (k,l)
    
    def plot_grid(self):
        for part in self.particles:
            color = [np.random.rand() for _ in range(3)]
            plt.fill_between([part.k-self.N/2, part.k-self.N/2+1], part.l-self.N/2, part.l-self.N/2+1, alpha=0.7, color=color) # complicated to explain this line, but just try it on paper to see it's correct
            plt.scatter(part.xh, part.yh, color='black')    
        edges = np.linspace(-self.N/2, self.N/2, self.N+1)
        plt.xticks(edges)
        plt.yticks(edges)
        plt.grid()
        plt.xlim([-self.N/2, self.N/2])
        plt.ylim([-self.N/2, self.N/2])
            
    # Given a particle calculates distance (in units of h) from lower left corner of cell c=(k,l) and the particle position
    # if c is not provided it uses the closest cell to the particle
    def get_eps(self, part, c=None):
        if c is None:
            c = [part.k, part.l]
        dist_x = part.xh-part.k+self.N/2
        dist_y = part.yh-part.l+self.N/2
        return (dist_x, dist_y)
    
    def top_hat(x):
        if abs(x) > 0.5:
            return 0
        else:
            return 1

    # calculates value of rho for cell (k,l)
    def calculate_rho(k,l):
        for part in self.particles:
            return 
    