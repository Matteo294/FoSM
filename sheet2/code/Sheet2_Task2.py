#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt


def Tf(T): # enter T in Kelvin!!; =^ rhs of ODE

    C = -4.83e-7 # Ks⁻¹
    if T <= 2e4:
        return C*(T/2e4)**10
    else: 
        return C*(T/2e4)**(-0.5) 


# 2nd order RK predictor corrector meth: K1 = f(y_n, t_n); K2 = f(y_n + K1*Delta_t, t_{n+1}); y_{n+1} = y_n + (K1 + K2)/2*Delta_t

def RK2pc(Deltat):

    #Deltat = 1e10 # s
    T = [1e7] # K
    t = [0] # s
    i = 0
    while T[i] >= 6e3:
        K1 = Tf(T[i])
        K2 = Tf(T[i] + K1*Deltat)
        Tnew = T[i] + (K1+K2)/2*Deltat
        T.append(Tnew)
        t.append(Deltat*(i+1))
        i += 1
    print("Predictor-corrector:", i, "steps")
    return T, t

# Adaptive stepszie RK2
def RK2as(Deltat, err_max):
    T = [1e7] # K
    t = [0] # s
    i = 0
    while T[i] >= 6e3:
        # Full Step
        K1 = Tf(T[i])
        K2 = Tf(T[i] + K1*Deltat)
        Tnew_full = T[i] + (K1+K2)/2*Deltat
        T.append(Tnew_full)
        t.append(t[-1]+Deltat)
        # Two half steps
        Deltat /= 2
        Tnew_half = T[i]
        for _ in range(2):
            K1 = Tf(T[i])
            K2 = Tf(T[i] + K1*Deltat)
            Tnew_half += (K1+K2)/2*Deltat
        err = abs(Tnew_full - Tnew_half)
        Deltat = 2*Deltat * (err_max/err)**(1/3)
        i += 1
    print("Adaptive stepsize:", i, "steps" )
    return T,t


if __name__ == "__main__":
    
    Deltats = [1e10, 2e10, 5e10, 8e10, 8.5e10, 8.8e10, 8.9e10, 9e10]

    for i in range(len(Deltats)):
        T, t = RK2pc(Deltats[i])
        j = 0
        #while T[j] > 2e4:
        #    j += 1
        fig = plt.figure()
        plt.plot(t,np.log(T),'-k', label = 'temperature evolution ln(T(t))')
        plt.xlabel('time (s)')
        plt.ylabel('ln(T)')
        plt.hlines(np.log(2e4),t[0],t[-1],colors='r',linestyles='--',label='$T=T_0$')
        #plt.plot(t[j],np.log(T[j]),'ro',label='$T=T_0$')
        plt.legend()
        plt.savefig(f'RK2pc_deltat{i}.png')
        #plt.show()   
        print(len(T), end=" ", flush=True)
    
    print("\n")
    
    for i in range(len(Deltats)):
        T, t = RK2as(Deltats[i], 50)
        j = 0
        #while T[j] > 2e4:
        #    j += 1
        fig = plt.figure()
        plt.plot(t,np.log(T),'-k', label = 'temperature evolution ln(T(t))')
        plt.xlabel('time (s)')
        plt.ylabel('ln(T)')
        plt.hlines(np.log(2e4),t[0],t[-1],colors='r',linestyles='--',label='$T=T_0$')
        #plt.plot(t[j],np.log(T[j]),'ro',label='$T=T_0$')
        plt.legend()
        plt.savefig(f'RK2as_deltat{i}_as.png')
        #plt.show()   
        print(len(T), end="", flush=True)
