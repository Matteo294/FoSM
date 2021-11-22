import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from Particle import *
import seaborn as sb

H = 4
N = 50
h = 2*H/N
nparticles = 1000

grid1 = Grid(H, N)
for n in range(nparticles):
    x, y = np.random.randn(2)
    if abs(x) < (H-h) and abs(y) < (H-h):
        grid1.add_particle([x,y])
grid1.plot_grid()
plt.show()

sb.heatmap(grid1.W0, cbar=False)
plt.savefig("images/zeroth_order.png", ppi=300)
sb.heatmap(grid1.W1, cbar=False)
plt.savefig("images/first_order.png", ppi=300)
sb.heatmap(grid1.W2, cbar=False)
plt.savefig("images/second_order.png", ppi=300)
plt.show()