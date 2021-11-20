import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from Particle import *
import seaborn as sb
############################### Exercise 1 #####################################
H = 4
N = 16
nparticles = 1000

grid1 = Grid(H, N)
for n in range(nparticles):
    grid1.add_particle(np.random.randn(10))
    print("Particle", n, "has coords", grid1.particles[-1].xh, grid1.particles[-1].yh)
    print("The closes cell is", grid1.particles[-1].k, grid1.particles[-1].l)
    print("The distance from the particle and the bottom left corner fo the closest cell is", grid1.get_eps(grid1.particles[-1]))
grid1.plot_grid()
plt.show()

sb.heatmap(grid1.W0)
plt.show()
################################################################################


############################### Exercise 2 #####################################
H = 4
N = 8
nparticles = 2

grid2 = Grid(H, N)
for n in range(nparticles):
    grid2.add_particle([(np.random.rand()*2-1)*H for _ in range(2)])
grid2.plot_grid()

################################################################################