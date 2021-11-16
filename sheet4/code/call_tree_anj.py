from tree_anj import *
from random import *
import time
from copy import *
from math import sqrt
import pandas as pd


# direct force calculation
def calc_direct(particles,nparticles):

    fexact = np.zeros((nparticles,3))
    for i in range(nparticles):
        xi = np.array(particles[i])
        for j in range(nparticles):
            if i == j:
                continue
            xj = np.array(particles[j])
            fexact[i,:] += -1/sqrt(np.linalg.norm(xi-xj)**2 + 0.001**2)**3*(xi-xj)
    return fexact


# generates N particles
def gen_parts(N):

    #
    # Create a set of randomly positioned particles
    # For convenience we assume all masses to be 1.
    # If we have in reality another mass, we can simply
    # rescale our answers.
    #
    nparticles = N
    particles = []
    for i in range(nparticles):
        x = random()
        y = random()
        z = random()
        particles.append([x,y,z])

    return particles


# tree based force calculation + returns runtime and average node interactions per particle
def calc_tree(particles,theta):
    
    #
    # Now create the tree
    #
    
    q=TreeClass(particles)
    q.insertallparticles()
    q.computemultipoles(0)

    t0 = time.time()
    q.allgforces(theta)
    t1 = time.time()
    runtime = t1-t0

    fapprox = deepcopy(q.forces)
    counts = np.array(q.counts)
    avc = np.mean(q.counts)

    return fapprox, runtime, avc


# calculates mean relative force deviation for given tree forces and direct sum forces
def calc_mean_eta(fapprox, fexact):
    diffs = fapprox - fexact
    diff_abs = np.linalg.norm(diffs,axis=1)

    eta = diff_abs/np.linalg.norm(fexact,axis=1)

    mean_eta = np.mean(eta)

    #print(eta)

    return mean_eta



# exercises b) and c)
def exbc():

    particles = gen_parts(100)

    print("starting tree gravity")

    fapprox, runtime, avc = calc_tree(particles,0.8)

    print("done in ", runtime, " seconds\n")
        
    print(f'average number of node interactions per particle: {avc}\n\n')

    #print(f'tree approximation force: {fapprox}')


    print("starting N^2 gravity")

    t0 = time.time()
    fexact = calc_direct(particles,100)
    t1 = time.time()
    fullgrav_dt = t1-t0
    #print(f'force from direct summation: {fexact}')
    print("done in ", fullgrav_dt, " seconds\n")

    #fexact = deepcopy(q.forces)

    # 
    # Now compare the approximate and exact versions

    fapprox = np.array(fapprox)
    mean_eta = calc_mean_eta(fapprox,fexact)
    print(f'<eta> = {mean_eta}')




# ex d) & e)

def exde():

    N = [5000,10000,20000,40000]
    #N = [5000]
    theta = [0.2, 0.4, 0.8]

    times_tree = np.zeros((len(N),len(theta)))
    times_direct = np.zeros(len(N))
    mean_etas = np.copy(times_tree)
    counts = np.copy(times_tree)

    for i in range(len(N)):
        print(f'nbr of particles: N = {N[i]}')
        particles = gen_parts(N[i])
        t0 = time.time()
        fexact = calc_direct(particles,N[i])
        t1 = time.time()
        runtime_direct = t1-t0
        times_direct[i] = runtime_direct
        for j in range(len(theta)):
            print(f'opening angle: theta = {theta[j]}')

            fapprox, runtime, avc = calc_tree(particles, theta[j])

            times_tree[i,j] = runtime
            counts[i,j] = avc

            fapprox = np.array(fapprox)

            mean_eta = calc_mean_eta(fapprox,fexact)

            mean_etas[i,j] = mean_eta

            print(f'\nruntimes: {runtime} vs. {runtime_direct}\n')
            print(f'\n<eta>: {mean_eta}\n')
            print(f'\ncounts: {avc}\n')


    dataframe1 = pd.DataFrame({'theta=0.2':times_tree[:,0],'theta=0.4':times_tree[:,1],'theta=0.8':times_tree[:,2]})
    dataframe1.to_csv('times_tree.dat')
    dataframe2 = pd.DataFrame({f'times for N = {N}':times_direct[:]})
    dataframe2.to_csv('times_direct.dat')
    dataframe3 = pd.DataFrame({'theta=0.2':mean_etas[:,0],'theta=0.4':mean_etas[:,1],'theta=0.8':mean_etas[:,2]})
    dataframe3.to_csv('mean_etas.dat')
    dataframe4 = pd.DataFrame({'theta=0.2':counts[:,0],'theta=0.4':counts[:,1],'theta=0.8':counts[:,2]})
    dataframe4.to_csv('counts.dat')
    print(f'runtimes tree: {times_tree}\n\n')
    print(f'runtimes direct: {times_direct}\n\n')
    print(f'<eta>: {mean_etas}\n\n')
    print(f'average node interactions: {counts}')
    
    

if __name__ == '__main__':

    exbc()
