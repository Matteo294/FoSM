import numpy as np
import matplotlib.pyplot as plt

def pnew(x):# normiert auf [1,20]
    return 20/19*1/x**2

def exact_inversion(x):
    return 1/(1 - 19/20*x)


def metropolis_mc(p,n,sigma=0.1):

    xstart = float(np.random.uniform(1,20,1))
    #xstart = 10
    x = []
    x.append(xstart)
    for i in range(n-1):
        p_i = p(xstart)
        x_ii = float(np.random.normal(x[-1],sigma))
        while x_ii < 1 or x_ii > 20 :
            x_ii = float(np.random.normal(x[-1],sigma))
        p_ii = p(x_ii)
        #print(p_ii/p_i)
        r = min(1,p_ii/p_i)
        print(r)
        u = np.random.uniform(0,1,1)
        if u <= r:
            x.append(x_ii)
        else:
            x.append(x[-1])

    return x

def ex3():

    x = metropolis_mc(pnew,int(1e6),sigma=0.1)

    #print(x)
    print(max(x))

    fig = plt.figure()
    bins = 100
    X = np.random.uniform(0,1,int(1e6))
    Y = exact_inversion(X)
    #plt.hist(Y,bins,density=True,color='r',alpha=0.5,label='exact inversion')
    plt.hist(x,bins,density=True,color='k',label='Metropolis Monte Carlo')
    plt.legend() 
    plt.show()


if __name__ == '__main__':
    ex3()

    
