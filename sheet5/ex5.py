import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from mpl_toolkits.mplot3d import Axes3D

h = 1
H = 4
N = 8
h = 2*H/N
nparticles = 3

# returns cell indices k,l for the particle lying at position r
def find_kl(r):
    if r[0] > 0:
        k = int(r[0]/h) + 1
    else:
        k = int(r[0]/h) - 1
    if r[1] > 0:
        l = int(r[1]/h) + 1
    else:
        l = int(r[1]/h) - 1
    return (k,l)

def plot(particles, H, N):
    h = 2*H/N
    for part in particles:
        color = [np.random.rand() for _ in range(3)]
        plt.fill_between([part.xk-0.5*h, part.xk+0.5*h], part.yl-0.5*h, part.yl+0.5*h, alpha=0.7, color=color)
        plt.scatter(part.x, part.y, color='black')    
    edges = np.linspace(-H, H, N+1)
    plt.xticks(edges)
    plt.yticks(edges)
    plt.grid()
    plt.xlim([-H, H])
    plt.ylim([-H, H])
    plt.show()


class Particle:

    def __init__(self, r):
        self.x = r[0]
        self.y = r[1]
        self.r = r
        self.k, self.l = find_kl(r)
        if(self.x >= 0):
            self.xk = self.k*h - 0.5
        else:
            self.xk = self.k*h + 0.5
        if(self.l >=0):
            self.yl = self.l*h - 0.5
        else:
            self.yl = self.l*h + 0.5
    
particles = []
for npart in range(nparticles):
    particles.append(Particle([(np.random.rand()*2-1)*H for _ in range(2)]))

plot(particles, H, N)