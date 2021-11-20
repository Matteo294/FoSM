import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from Particle import *
import seaborn as sb

H = 4
N = 32
h = 2*H/N
nparticles = 1000

grid1 = Grid(H, N)
for n in range(nparticles):
    x, y = np.random.randn(2)
    if abs(x) < (H-1.1*h) and abs(y) < (H-1.1*h):
        grid1.add_particle([x,y])
        print("Particle", n, "has coords", grid1.particles[-1].xh, grid1.particles[-1].yh)
        print("The closes cell is", grid1.particles[-1].k, grid1.particles[-1].l)
        print("The distance from the particle and the bottom left corner fo the closest cell is", grid1.get_eps(grid1.particles[-1]))
grid1.plot_grid()
plt.show()

plt.subplot(1,2,1)
sb.heatmap(grid1.W0)
plt.subplot(1,2,2)
sb.heatmap(grid1.W1)
plt.show()