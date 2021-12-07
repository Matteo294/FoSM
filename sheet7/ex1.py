import numpy as np
import matplotlib.pyplot as plt

# a)
def randu_gen(I=1,n=1):
    """I : seed, needs to be an odd integer; n : number of recursions, i.e. number of output random numbers // returns lists I_n, u_n lists of generated random numbers as integers I_n and normalized to [0,1] u_n
    """
    # sth like I = 2*I + 1 to ensure it is an odd integer?
    I_n = [I]
    u_n = []
    for i in range(n):
        II = np.double((65539*I_n[-1])%2e31) # where should I use double?
        uii = II/2.e31
        I_n.append(II)
        u_n.append(uii)
    return I_n[1:], u_n # do I want to include or exclude seed?

def exb():
    
    print('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    print('Exercise 1.b),c):')
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    N = int(input('number of points N = '))
    I = int(input('seed = ')) 
    print(f'\nstarting seed = {I}')
    tuples = np.zeros((N,2))
    # or instead create randu(I,n=2000) and use every two successive pairs as points? results in the same
    # or use list of 1000 seeds for every pair?
    #for i in range(N):
    #    I_2, u_2 = randu_gen(I,n=2)
    #    tuples[i,:] = u_2
    #    I = 2*I_2[-1] + 1
    I_2, u_2 = randu_gen(I,n=2*N)
    tuples[:,0] = u_2[0::2]
    tuples[:,1] = u_2[1::2]


    fig0 = plt.figure()
    ax0 = fig0.add_subplot()
    ax0.scatter(tuples[:,0],tuples[:,1],c='k',marker='.',label='successive random numbers')
    ax0.set_xlabel(r'$x_i\equiv u_{2i}$')
    ax0.set_ylabel(r'$y_i\equiv u_{2i+1}$')
    ax0.set_title('Tuples of successive random numbers')

    triples = np.zeros((N,3))
    #for i in range(N):
    #    I_3, u_3 = randu_gen(I,n=3)
    #    triples[i,:] = u_3
    #    I = 2*I_3[-1] + 1
    I_3, u_3 = randu_gen(I,n=3*N)
    for i in range(3):
        triples[:,i] = u_3[i::3]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(projection='3d')
    ax1.scatter(triples[:,0],triples[:,1],triples[:,2],c='k',marker='.',label='successive random numbers')
    ax1.set_xlabel(r'$x_i\equiv u_{3i}$')
    ax1.set_ylabel(r'$y_i\equiv u_{3i+1}$')
    ax1.set_zlabel('$z_i\equiv u_{3i+2}$')
    ax1.set_title('Triples of successive random numbers')
    
    plt.show()

def exc():

    print('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    print('Exercise 1.c):')
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    
    xmin = 0.2
    xmax = 0.201
    ymin = 0.3
    ymax = 0.301

    X = []
    Y = []

    I = 111
    print(f'\nstarting seed = {I}')
    
    while len(X) < 1000:        
        I_2, u_2 = randu_gen(I,n=2) 
        if xmin <= u_2[0] <= xmax and ymin <= u_2[1] <= ymax:
            X.append(u_2[0])
            Y.append(u_2[1])
        I = 2*I_2[-1] + 1
        print(len(X)) 
    fig0 = plt.figure()
    ax0 = fig0.add_subplot()
    ax0.scatter(X,Y,c='k',marker='.',label='successive random numbers')
    ax0.set_xlabel(r'$x_i\equiv u_{2i}$')
    ax0.set_ylabel(r'$y_i\equiv u_{2i+1}$')
    ax0.set_title('Tuples of successive random numbers')
    plt.show()

def exd():

    
    print('\n ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    print('Exercise 1.d):')
    print('\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')

    N = int(input('number of points N = '))
    #I = int(input('seed = ')) 
    #print(f'\nstarting seed = {I}')
    tuples = np.zeros((N,2))
    u_2 = np.random.uniform(size=2*N) 
    tuples[:,0] = u_2[0::2]
    tuples[:,1] = u_2[1::2]


    fig0 = plt.figure()
    ax0 = fig0.add_subplot()
    ax0.scatter(tuples[:,0],tuples[:,1],c='k',marker='.',label='successive random numbers')
    ax0.set_xlabel(r'$x_i\equiv u_{2i}$')
    ax0.set_ylabel(r'$y_i\equiv u_{2i+1}$')
    ax0.set_title('Tuples of successive random numbers using np.random.uniform')

    triples = np.zeros((N,3))
    u_3 = np.random.uniform(size=3*N) 
    for i in range(3):
        triples[:,i] = u_3[i::3]

    fig1 = plt.figure()
    ax1 = fig1.add_subplot(projection='3d')
    ax1.scatter(triples[:,0],triples[:,1],triples[:,2],c='k',marker='.',label='successive random numbers')
    ax1.set_xlabel(r'$x_i\equiv u_{3i}$')
    ax1.set_ylabel(r'$y_i\equiv u_{3i+1}$')
    ax1.set_zlabel('$z_i\equiv u_{3i+2}$')
    ax1.set_title('Triples of successive random numbers')
    
    plt.show()

if __name__ == "__main__":
    #exb()
    exd()
    #exc()
