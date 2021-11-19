import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from Particle import *

############################### Exercise 1 #####################################
H = 4
N = 8
nparticles = 20

grid1 = Grid(H, N)
for n in range(nparticles):
    grid1.add_particle([(np.random.rand()*2-1)*H for _ in range(2)])
grid1.plot_grid()
################################################################################