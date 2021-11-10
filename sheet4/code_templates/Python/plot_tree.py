import numpy as np
import pylab as plt
from mpl_toolkits.mplot3d import Axes3D
import os as os
from tree import *
from random import *

#======================================================================
#                            PLOT TREE
#                         by Paul Molliere 
#                              2014
#         (this version uses tree.py by Cornelis Dullemond)
#======================================================================

plt.ion()

fig = plt.figure()
ax = fig.gca(projection='3d')

def plot_cube(x,y,z,size):
    x12 = np.array([x-size/2.,x+size/2.])
    y12 = np.array([y-size/2.,y+size/2.])
    z12 = np.array([z-size/2.,z+size/2.])
    ones = np.ones(2)
    ax.plot(x12,y12[0]*ones,z12[0]*ones,'-',color='red')
    ax.plot(x12,y12[1]*ones,z12[0]*ones,'-',color='red')
    ax.plot(x12,y12[0]*ones,z12[1]*ones,'-',color='red')
    ax.plot(x12,y12[1]*ones,z12[1]*ones,'-',color='red')
    ax.plot(x12[0]*ones,y12,z12[0]*ones,'-',color='red')
    ax.plot(x12[1]*ones,y12,z12[0]*ones,'-',color='red')
    ax.plot(x12[0]*ones,y12,z12[1]*ones,'-',color='red')
    ax.plot(x12[1]*ones,y12,z12[1]*ones,'-',color='red')
    ax.plot(x12[0]*ones,y12[0]*ones,z12,'-',color='red')
    ax.plot(x12[0]*ones,y12[1]*ones,z12,'-',color='red')
    ax.plot(x12[1]*ones,y12[0]*ones,z12,'-',color='red')
    ax.plot(x12[1]*ones,y12[1]*ones,z12,'-',color='red')

nparticles = 30
particles = []
for i in range(nparticles):
    x = random()
    y = random()
    z = random()
    particles.append([x,y,z,1.0])

q=TreeClass(particles)
q.insertallparticles()
q.computemultipoles(0)

masses = np.array(q.particles)
plt.plot(masses[:,0],masses[:,1],masses[:,2],'o',color='blue')

for i in range(len(q.nodelist)):
    plot_cube(q.nodelist[i].xc[0],q.nodelist[i].xc[1],q.nodelist[i].xc[2],q.nodelist[i].size)

# UNCOMMENT THE FOLLOWING BLOCK TO CREATE A MOVIE OF YOUR OCT-TREE (TAKES SOME TIME) 
# for ii in range(0,360,1):
#    print(ii)
#    ax.view_init(elev=10., azim=ii)
#    plt.savefig("tree_plots/movie"+str(ii).zfill(3)+".png")
# os.system('ffmpeg -i tree_plots/movie%03d.png tree_code.mov')

