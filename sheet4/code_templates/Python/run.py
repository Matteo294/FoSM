from tree import *
from random import *
import numpy as np
from matplotlib import pyplot as plt

# computes force felt by x because of y
def compute_force(x, y, m1=1, m2=1):
    r2 = np.linalg.norm((x-y))
    return -m1*m2*(x-y)/r2**3

nparticles = 100
particles = []
for i in range(nparticles):
    x = random()
    y = random()
    z = random()
    particles.append(np.array([x,y,z]))
q=TreeClass(particles)
q.insertallparticles()
q.computemultipoles(0)
q.allgforces(0.8)

err = np.zeros((nparticles, 3))
force = np.zeros((nparticles, 3))
for i in range(nparticles):
    force[i][:] = sum([compute_force(particles[i], particles[j]) for j in range(nparticles) if i != j])
    err[i][:] = abs((q.forces[i][:] - force[i][:])/force[i][:])

nrange = range(nparticles)
plt.plot(nrange, err[:,0])
plt.plot(nrange, err[:,1])
plt.plot(nrange, err[:,2])
plt.show()