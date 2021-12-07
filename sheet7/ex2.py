import numpy as np
import matplotlib.pyplot as plt
from time import process_time

def func(x):
    
    d = np.shape(x)[-1]
    f = 1
    for i in range(d):
        f *= 3/2*(1-x[i]**2)
    return f

def midpoint(f,d,n):

    m = np.arange(1/(2*n),1+1/(2*n),1/n)
    M = np.array(np.meshgrid(*(m,)*d)).T.reshape(-1,d)
    D = np.shape(M)[0]
    fM = np.zeros(D)
    for i in range(D):
        fM[i] = f(M[i])
    # I = sum_i^n**2 1/n**2 f(m_i) where m_i is midpoint in each subcube
    I = 1/n**2*np.sum(fM)
    return I

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

    for i in d:
        start = process_time()
        Imp = midpoint(func,i,n)
        stop = process_time()
        t_mpoint.append(stop-start)
        I_mpoint.append(Imp)

        start = process_time()
        Imc = monte_carlo_int(func,i,N)
        stop = process_time()
        t_mc.append(stop-start)
        I_mc.append(Imc)

    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('\n Exercise 2\n')
    print('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~')
    print('\n Midpoint method:\n')
    print(f'I for d in {1,2,...,10}: \n')
    print(f'I = {I_mpoint}')
    print('\n CPU times for each dimension: \n')
    print(f't = {t_mpoint}')

    print('\n Monte Carlo method:\n')
    print(f'I for d in {1,2,...,10}: \n')
    print(f'I = {I_mc}')
    print('\n CPU times for each dimension: \n')
    print(f't = {t_mc}')

if __name__ == '__main__':
    ex2()



        
