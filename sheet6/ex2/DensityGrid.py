from matplotlib import pyplot as plt 
import numpy as np
import seaborn as sb

############################# Grid class ###########################
class DensityGrid:
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
        def __init__(self, r, grid_instance, mass=1):
            selfdensity_grid = grid_instance
            self.x = r[0] # x coord
            self.y = r[1] # y coord
            self.r = r # position vector
            self.k, self.l = selfdensity_grid.find_kl(self.r) # get indices of the closest cell
            self.xh = self.x/selfdensity_grid.h # x coordinate in units of h
            self.yh = self.y/selfdensity_grid.h # y coordinate in units of h
            self.rg = [self.xh, self.yh] # r in units of h
            self.m = mass

    # Adds a particle to the grid at position r
    def add_particle(self, r):
        self.particles.append(self.Particle(r, self))
        k = self.particles[-1].k
        l = self.particles[-1].l
        m = self.particles[-1].m
        # Update weights for zeroth order
        self.W0[k,l] += m
        # Update weights for the first order
        epsx, epsy = self.get_eps(self.particles[-1])
        self.W1[k,l] += m*self.w1_func(epsx)*self.w1_func(epsy)
        self.W1[k,l+1] += m*self.w1_func(epsx)*self.w1_func(3/2-epsy)
        self.W1[k,l-1] += m*self.w1_func(epsx)*self.w1_func(1/2+epsy)
        self.W1[k+1,l] += m*self.w1_func(3/2-epsx)*self.w1_func(epsy)
        self.W1[k+1,l+1] += m*self.w1_func(3/2-epsx)*self.w1_func(3/2-epsy)
        self.W1[k+1,l-1] += m*self.w1_func(3/2-epsx)*self.w1_func(1/2+epsy)
        self.W1[k-1,l] += m*self.w1_func(1/2+epsx)*self.w1_func(epsy)
        self.W1[k-1,l+1] += m*self.w1_func(1/2+epsx)*self.w1_func(3/2-epsy)
        self.W1[k-1,l-1] += m*self.w1_func(1/2+epsx)*self.w1_func(1/2+epsy)
        # Update weights for the second order
        d1x = epsx - 1/2 # distance from neirest cell
        d2x = 3/2 - epsx # distance from next cell
        d0x = 1/2 + epsx # distance from previous cell
        d1y = epsy - 1/2 
        d2y = 3/2 - epsy
        d0y = 1/2 + epsy
        self.W2[k,l] += m*self.w2_func(d1x)*self.w2_func(d1y)
        self.W2[k,l+1] += m*self.w2_func(d1x)*self.w2_func(d2y)
        self.W2[k,l-1] += m*self.w2_func(d1x)*self.w2_func(d0y)
        self.W2[k+1,l] += m*self.w2_func(d2x)*self.w2_func(d1y)
        self.W2[k+1,l+1] += m*self.w2_func(d2x)*self.w2_func(d2y)
        self.W2[k+1,l-1] += m*self.w2_func(d2x)*self.w2_func(d0y)
        self.W2[k-1,l] += m*self.w2_func(d0x)*self.w2_func(d1y)
        self.W2[k-1,l+1] += m*self.w2_func(d0x)*self.w2_func(d2y)
        self.W2[k-1,l-1] += m*self.w2_func(d0x)*self.w2_func(d0y)

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
        if len(self.particles) < 20:
            plt.xticks(edges, labels=edges)
            plt.yticks(edges, labels=edges)
        plt.grid()
        plt.xlim([-self.N/2, self.N/2])
        plt.ylim([-self.N/2, self.N/2])
        plt.xlabel("x/h")
        plt.ylabel("y/h")
            
    # Given a particle calculates distance (in units of h) from lower left corner of cell c=(k,l) and the particle position
    # if c is not provided it uses the closest cell to the particle
    def get_eps(self, part, c=None):
        if c is None:
            c = [part.k, part.l]
        dist_x = part.xh-c[0]+self.N/2
        dist_y = part.yh-c[1]+self.N/2
        return float(dist_x), float(dist_y)
    
    def top_hat(self, x):
        if abs(x) > 0.5:
            return 0
        else:
            return 1
    
    def w1_func(self, x):
        if abs(x) < 1/2:
            return 1/2*(1 - abs(x))
        else:
            return 0

    def w2_func(self, x):
        if abs(x) > 3/2:
            return 0
        elif abs(x) < 3/2 and abs(x) > 1/2:
            return 1/2*(3/2-abs(x))**2
        else:
            return (3/4-x**2)
    