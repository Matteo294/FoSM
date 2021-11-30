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
    x, y = np.random.uniform(low=-grid.H + grid.h, high = grid.H - grid.h, size=2)
    X = np.array((x,y))

    print(X)

    grid.add_particle(X)

    rho_real = grid.rho1
    

    # b)

    #kx = 2*np.pi/L*np.arange(-N/2,N/2,1) # !! N must be even for this to work
    kx = 2*np.pi/L*np.arange(0,N,1)
    ky = np.flip(kx)
    #k = np.vstack((kx,ky))
    #print(kx)
    
    greens_kspace = np.zeros((len(kx),len(ky)))
    for i in range(len(kx)):
        for j in range(len(ky)):
            k = np.array((kx[i],ky[j]))
            greens_kspace[j,i] = greens_func_kspace(k) 

    #print(greens_kspace)
    greens_kspace_padded = np.zeros((2*N,2*N))
    greens_kspace_padded[0:N,0:N] = greens_kspace_padded[0:N,N:2*N] = greens_kspace_padded[N:2*N,0:N] = greens_kspace_padded[N:2*N,N:2*N] = greens_kspace

    #greens_kspace_padded = np.fft.fftshift(greens_kspace_padded)
    rho_real_padded = np.zeros((2*N,2*N))
    rho_real_padded[0:N,0:N] = rho_real # rho lives in lower left corner
    rho_kspace_padded = np.fft.fft2(rho_real_padded)
    
    phi_kspace_padded = rho_kspace_padded*greens_kspace_padded

    phi_real_padded = np.fft.ifft2(phi_kspace_padded).real
    #print(phi_real_padded)
    grid.show_density_map(grid.rho1)
    phi_real = phi_real_padded[0:N,0:N]
    #print(phi_real)

    # c)
   
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

    fig0 = plt.figure()
    Z0 = plt.imshow(A[0])
    fig0.colorbar(Z0)
    fig1 = plt.figure()
    Z1 = plt.imshow(A[1]) 
    fig1.colorbar(Z1)   
    plt.show()


    # d)
    
    rmax = L/2
    rmin = 0.3*L/N
    rrel = rmax/rmin

    positions = []
    dists = []
    for i in range(100):
        p, q = np.random.uniform(0,1,size=2)
        dx = rmin*rrel**p*np.cos(2*np.pi*q)
        dy = rmin*rrel**p*np.sin(2*np.pi*q)
        r = np.sqrt(dx**2 + dy**2)
        if -L/2+grid.h <= X[0] + dx < L/2-grid.h: 
            x = X[0] + dx
        else:
            x = X[0] - dx
        if -L/2+grid.h <= X[1] + dy < L/2-grid.h:
            y = X[1] + dy
        else:
            y = X[1] - dy
        pos = np.array((x,y))    
        dists.append(r)
        positions.append(pos)

    positions = np.array(positions)
    fig = plt.figure()
    plt.plot(X[0],X[1],'xr')
    plt.plot(positions[:,0],positions[:,1],'bo')
    plt.xlim(-0.5,0.5)
    plt.ylim(-0.5,0.5)
    plt.show()

    # intepolate force by CIC
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
    abs_accs = np.linalg.norm(accs, axis=-1)
    
    fig = plt.figure()
    plt.plot(dists,abs_accs,'k.')
    xx = np.linspace(0.01,0.5,1000)
    plt.plot(xx,2/xx,'r-')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()

if __name__ == "__main__":
    ex2()





#TODOS:
# solve boundary issues
# check right fft / shifts etc
# check offset in a power law?
# search for a way that we can allow positions also in boundary cells --> resolve issues in W calculation
