import numpy as np
import matplotlib.pyplot as plt

def func(x):
    
    d = len(x)
    f = 1
    for i in range(d):
        f *= 3/2*(1-x[i]**2)
    return f

def midpoint(f,d,n):

    # I = sum_i^n**2 1/n**2 f(m_i) where m_i is midpoint in each subcube
    # midpoint integral, n subdivisions, d dimensions
    # create lists of x and y --> midpoints, evaluate f at each midpoint and sum up, np.sum matrix?
    pass

def monte_carlo_int(f,d,N):

    NN = np.random.uniform(size=N*d)
    X = NN.reshape((N,d))
    # monte carlo: I = V/N sum_i^N f(x_i), generate N vectors X of dimension d
    I = 0
    for i in range(N):
        I += f(X[i,:])
    I /= N
    return I

def ex2():
    
    n = 6
    N = 20000
    d = np.arange(1,10,1)
    t_mpoint = []
    t_mc = []
    I_mpoint = []
    I_mc = []

    pass
