import numpy as np
from Density_map import *

# 2)

L = 1
N = 256 

def greens_func_kspace(k):

    if k.all() == 0:
        return 0
    else:
        return -4*np.pi/(k[0]**2 + k[1]**2)



def ex2():

    # a)

    grid = Density_map(L/2,N)
    x, y = np.random.uniform(low=-grid.H + grid.h, high = grid.H - grid.h, size=2) # particle inside one cell length, aka not in a boundary cell
    X = np.array((x,y))

    print(f'particle at cartesian position: X = {X}')

    grid.add_particle(X)

    rho_real = grid.rho1 # denisty mesh as implemented in last weeks exercise
    

    # b)

    #kx = 2*np.pi/L*np.arange(-N/2,N/2,1) # !! N must be even for this to work
    kx = 2*np.pi/L*np.arange(0,N,1) # which one is correct? both should work equally? but the symmetric one gives weird chess board patterns in force?
    ky = np.flip(kx) # bc array and physical directions in y opposite
    #print(kx)
    
    greens_kspace = np.zeros((len(kx),len(ky))) # set up greens function at each cell position
    for i in range(len(kx)):
        for j in range(len(ky)):
            k = np.array((kx[i],ky[j]))
            greens_kspace[j,i] = greens_func_kspace(k) # x and y indices switched!!

    #print(greens_kspace)
    greens_kspace_padded = np.zeros((2*N,2*N)) # zero padding to prevent from feeling force from neighbouring periodicly continued cell
    greens_kspace_padded[0:N,0:N] = greens_kspace_padded[0:N,N:2*N] = greens_kspace_padded[N:2*N,0:N] = greens_kspace_padded[N:2*N,N:2*N] = greens_kspace # greens should satisfy periodicity in each quarter

    #greens_kspace_padded = np.fft.fftshift(greens_kspace_padded) # necessary or not?
    rho_real_padded = np.zeros((2*N,2*N))
    rho_real_padded[0:N,0:N] = rho_real # rho lives in lower left corner
    rho_kspace_padded = np.fft.fft2(rho_real_padded) # DFT of zero padded rho
    
    phi_kspace_padded = rho_kspace_padded*greens_kspace_padded 

    phi_real_padded = np.fft.ifft2(phi_kspace_padded).real # should yield right result in lower left corner for phi, discard other quarters
    #print(phi_real_padded)

    #grid.show_density_map(grid.rho1)

    phi_real = phi_real_padded[0:N,0:N] # lower left corner
    #print(phi_real)


    # c)

    # calculate force in each cell   
    A = np.zeros((2,N,N))
    for i in range(1,N-1):
        # central difference in x direc
        A[0,:,i] = -(phi_real[:,i+1]-phi_real[:,i-1])/(2*grid.h)
        # central diff in y direc
        A[1,i,:] = -(phi_real[i-1,:] - phi_real[i+1,:])/(2*grid.h)

    # at the boundaries:

    # use single step, i.e. BE and FE
    A[0,:,0] = -(phi_real[:,1]-phi_real[:,0])/grid.h 
    A[0,:,-1] = -(phi_real[:,-1]-phi_real[:,-2])/grid.h
    A[1,0,:] = -(phi_real[0,:]-phi_real[1,:])/grid.h
    A[1,-1,:] = -(phi_real[-2,:]-phi_real[-1,:])/grid.h

    #print(A)
    # visualize A mesh x and y components
    fig0 = plt.figure()
    Z0 = plt.imshow(A[0])
    fig0.colorbar(Z0)
    plt.title('$A_x$')
    fig1 = plt.figure()
    Z1 = plt.imshow(A[1]) 
    fig1.colorbar(Z1)
    plt.title('$A_y$')   
    plt.show()


    # d)
    
    rmax = L/2
    rmin = 0.3*L/N
    rrel = rmax/rmin
    # set up 100 randomly positioned evaluation points in given region around particle from a)
    positions = [] # stores positions
    dists = [] # stores absolute distance from particle
    for i in range(100):
        p, q = np.random.uniform(0,1,size=2)
        dx = rmin*rrel**p*np.cos(2*np.pi*q)
        dy = rmin*rrel**p*np.sin(2*np.pi*q)
        r = np.sqrt(dx**2 + dy**2)
        if -L/2+grid.h <= X[0] + dx <= L/2-grid.h: # for that we can use our method to setup W matrix for the positions, make sure that it's inside at least one cell length from boundary
            x = X[0] + dx
        else:
            x = X[0] - dx
        if -L/2+grid.h <= X[1] + dy <= L/2-grid.h:
            y = X[1] + dy
        else:
            y = X[1] - dy
        pos = np.array((x,y))    
        dists.append(r)
        positions.append(pos)

    positions = np.array(positions)

    # shows positions relative to particle
    fig = plt.figure()
    plt.plot(positions[:,0],positions[:,1],'bo')
    plt.plot(X[0],X[1],'xr')
    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.show()

    # intepolate force by CIC --> for every positions, force contributions from neighbouring cells by assignement function
    accs = []
    for i in range(100):
        grid_temporary = Density_map(L/2,N)
        grid_temporary.add_particle(positions[i])
        W = grid_temporary.W1
        #print(W)
        ax = np.sum(A[0,:,:]*W)
        ay = np.sum(A[1,:,:]*W)
        a = np.array((ax,ay))
        accs.append(a)
        #print(a)

    dists = np.array(dists)
    accs = np.array(accs)
    #print(np.shape(dists))
    #print(np.shape(accs))
    abs_accs = np.linalg.norm(accs, axis=-1) # (a_x^2 + a_y^2)^1/2 for each position 
    
    # plot absolute a as function of distance r from particle / comment out for e) when running 10 cycles
    #fig = plt.figure()
    #plt.plot(dists,abs_accs,'k.',label='experimental a values')
    #xx = np.linspace(0.3*L/N,0.5,1000)
    #plt.plot(xx,2/xx,'r-',label=r'power law $a \propto \frac{2}{r}$') # power law a \propto 2/r
    #plt.xscale('log')
    #plt.xlabel('r')
    #plt.yscale('log')
    #plt.ylabel('a')
    #plt.vlines(L/N, np.min(abs_accs)-0.3,2/xx[0]+10.3,color='g',linestyles='--',label='L/N')
    #plt.legend(loc='lower left')
    #plt.show()

    return dists, abs_accs 


def ex_e():
    dists_all = []
    accs_abs_all = []
    for i in range(10):
        dists, abs_accs = ex2()
        dists_all.append(dists)
        accs_abs_all.append(abs_accs)
    
    fig = plt.figure()
    plt.plot(dists_all[0],accs_abs_all[0],'k.',label='experimental a values')
    plt.plot(dists_all[1:],accs_abs_all[1:],'k.')
    xx = np.linspace(0.3*L/N,0.5,1000)
    plt.plot(xx,2/xx,'r-', label=r'power law $a\propto \frac{2}{r} $')
    plt.xscale('log')
    plt.xlabel('r')
    plt.yscale('log')
    plt.ylabel('a')
    plt.vlines(L/N, np.min(accs_abs_all)-0.3,2/xx[0]+10.3,color='g',linestyles='--',label='L/N')
    plt.legend(loc='lower left')
    plt.show()

if __name__ == "__main__":
    ex_e()
    #ex2()




#TODOS:
# solve boundary issues
# check right fft / shifts etc
# check offset in a power law?
# search for a way that we can allow positions also in boundary cells --> resolve issues in W calculation
