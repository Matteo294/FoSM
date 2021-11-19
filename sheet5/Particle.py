from matplotlib import pyplot as plt 
import numpy as np

###################### Grid class ########################
# Creates a 2D grid from -H to H divided into N^2 cells
# The side length of one square is h
# add particles to the grid with the subclass Particle
##########################################################
class Grid:
    def __init__(self, H, N):
        self.H = H 
        self.N = N
        self.h = 2*H/N
        self.particles = []

    # Adds a particle to the grid at position r
    def add_particle(self, r):
        self.particles.append(self.Particle(r, self))

    # Get indices of the cell closest to position r
    def find_kl(self, r):
        if r[0] > 0:
            k = int(r[0]/self.h) + 1
        else:
            k = int(r[0]/self.h) - 1
        if r[1] > 0:
            l = int(r[1]/self.h) + 1
        else:
            l = int(r[1]/self.h) - 1
        return (k,l)
    
    def plot_grid(self):
        for part in self.particles:
            color = [np.random.rand() for _ in range(3)]
            plt.fill_between([part.k, part.k-part.sgn_x], part.l, part.l-part.sgn_y, alpha=0.7, color=color) # complicated to explain this line, but just try it on paper to see it's correct
            plt.scatter(part.xh, part.yh, color='black')    
        edges = np.linspace(-self.H, self.H, self.N+1)
        plt.xticks(edges)
        plt.yticks(edges)
        plt.grid()
        plt.xlim([-self.H, self.H])
        plt.ylim([-self.H, self.H])
        plt.show()

    class Particle:
        def __init__(self, r, grid_instance):
            self.grid = grid_instance
            self.x = r[0] # x coord
            self.y = r[1] # y coord
            self.r = r # position vector
            self.k, self.l = self.grid.find_kl(self.r) # get indices of the closest cell
            self.xh = self.x/self.grid.h # x coordinate in units of h
            self.yh = self.y/self.grid.h # y coordinate in units of h
            # The following variables keep track of the sign of the coordinates (will be useful)
            self.sgn_x = 1
            self.sgn_y = 1
            if self.x < 0:
                self.sgn_x = -1
            if self.y < 0:
                self.sgn_y = -1
            