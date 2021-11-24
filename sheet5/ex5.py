import numpy as np
from matplotlib import pyplot as plt
from numpy.random.mtrand import random
from DensityGrid import *
import seaborn as sb

############################## Exercises 2.1 and 2.2 ####################################
H = 10
N = 10
h = 2*H/N
nparticles = 1
grid = DensityGrid(H, N)
for n in range(nparticles):
    x, y = np.random.randn(2)
    while abs(x) > (H-h) or abs(y) > (H-h):
        x, y = np.random.randn(2)
    grid.add_particle([x,y])
    print("Particle", n, "is at:", grid.particles[-1].xh, grid.particles[-1].yh)
    print("The closest cell is:", grid.particles[-1].k, grid.particles[-1].l)

grid.plot_grid()
plt.show()


# Zeroth order
ax = sb.heatmap(grid.W0.T, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.title("N = 1")
plt.savefig("./images/1p_zeroth_order.eps", dpi=300)
plt.show()

# First order
ax = sb.heatmap(grid.W0.T, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.title("N = 1")
plt.savefig("./images/1p_first_order.eps", dpi=300)
plt.show()

# Second order
ax = sb.heatmap(grid.W2.T, vmin=0.0,cmap='inferno')
ax.invert_yaxis()
plt.title("N = 1")
plt.savefig("./images/1p_second_order.eps", dpi=300)
plt.show()

print("M zeroth order:", np.sum(np.sum(grid.W0)))
print("M first order: ", np.sum(np.sum(grid.W1)))
print("M second order: ", np.sum(np.sum(grid.W2)))
print("\n")




############################## Exercise 2.3 ####################################
H = 10
N = 40
h = 2*H/N
nparticles = 100
grid = DensityGrid(H, N)
for n in range(nparticles):
    x, y = np.random.normal(loc=0.0, scale=3.0, size=2)
    if abs(x) < (H-h) and abs(y) < (H-h):
        x, y = np.random.randn(2)
        grid.add_particle([x,y])

# Zeroth order
ax = sb.heatmap(grid.W0, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.xlim([-5/h+N/2, 5/h+N/2])
plt.ylim([-5/h+N/2, 5/h+N/2])
plt.title("N = 100")
plt.savefig("./images/100p_zeroth_order.eps", dpi=300)
plt.show()

# First order
ax = sb.heatmap(grid.W1, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.xlim([-5/h+N/2, 5/h+N/2])
plt.ylim([-5/h+N/2, 5/h+N/2])
plt.title("N = 100")
plt.savefig("./images/100p_first_order.eps", dpi=300)
plt.show()

# Second order
ax = sb.heatmap(grid.W2, vmin=0.0,cmap='inferno')
ax.invert_yaxis()
plt.xlim([-5/h+N/2, 5/h+N/2])
plt.ylim([-5/h+N/2, 5/h+N/2])
plt.title("N = 100")
plt.savefig("./images/100p_second_order.eps", dpi=300)
plt.show()

print("M zeroth order:", np.sum(np.sum(grid.W0)))
print("M first order: ", np.sum(np.sum(grid.W1)))
print("M second order: ", np.sum(np.sum(grid.W2)))
print("Total number of particles:", len(grid.particles))
print("\n")


############################## Exercise 2.4 ####################################
H = 10
N = 100
h = 2*H/N
nparticles = 10000
grid = DensityGrid(H, N)
for n in range(nparticles):
    x, y = np.random.normal(loc=0.0, scale=3.0, size=2)
    if abs(x) < (H-h) and abs(y) < (H-h):
        x, y = np.random.randn(2)
        grid.add_particle([x,y])

# Zeroth order
ax = sb.heatmap(grid.W0, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.xlim([-4/h+N/2, 4/h+N/2])
plt.ylim([-4/h+N/2, 4/h+N/2])
plt.savefig("./images/zeroth_order.eps", dpi=300)
plt.title("N = 10000")
plt.show()

# First order
ax = sb.heatmap(grid.W1, vmin=0.0, cmap='inferno')
ax.invert_yaxis()
plt.xlim([-4/h+N/2, 4/h+N/2])
plt.ylim([-4/h+N/2, 4/h+N/2])
plt.savefig("./images/first_order.eps", dpi=300)
plt.title("N = 10000")
plt.show()

# Second order
ax = sb.heatmap(grid.W2, vmin=0.0,cmap='inferno')
ax.invert_yaxis()
plt.xlim([-4/h+N/2, 4/h+N/2])
plt.ylim([-4/h+N/2, 4/h+N/2])
plt.savefig("./images/second_order.eps", dpi=300)
plt.title("N = 10000")
plt.show()

print("M zeroth order:", np.sum(np.sum(grid.W0)))
print("M first order: ", np.sum(np.sum(grid.W1)))
print("M second order: ", np.sum(np.sum(grid.W2)))
print("Total number of particles:", len(grid.particles))
print("\n")





# weights functions plots
def zeroth_order_weightf(x):
    if abs(x) > 1/2:
        return 0
    else:
        return 1

def first_order_weightf(x):
    if abs(x) > 1:
        return 0
    else:
        return (1-abs(x))

xvals = np.linspace(-2, 2, 1000)
plt.plot(xvals, [zeroth_order_weightf(x) for x in xvals], label='zeroth', lw=2, c='black')
plt.plot(xvals, [first_order_weightf(x) for x in xvals], label='first', lw=2, c='red')
plt.plot(xvals, [grid.w2_func(x) for x in xvals], label='second', lw=2, c='blue')
plt.xlabel(r'$s=|x_i-x_p|$', fontsize=14)
plt.ylabel(r'$w(s)$', fontsize=14)
plt.legend()
plt.savefig("./images/weight_functions.pdf")
plt.show()