import numpy as np
from matplotlib import pyplot as plt
from DensityGrid import *

H = 1
N = 16
h = 2*H/N
nparticles = 1000
grid = DensityGrid(H, N)
for n in range(nparticles):
    x, y = [np.random.rand()*2 - 1 for _ in range(2)]
    while abs(x) > (H-h) or abs(y) > (H-h):
        x, y = [np.random.rand()*2 - 1 for _ in range(2)]
    grid.add_particle([x,y])

# First order
ax = sb.heatmap(grid.W1.T, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.show()