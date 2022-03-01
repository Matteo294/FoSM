import numpy as np 
from matplotlib import pyplot as plt

data = np.load("output.npz")

radii = [np.sqrt(x[0]**2+ x[1]**2 + x[2]**2) for x in data['x']]

rho = np.zeros(data['x'].shape[0])
for i in range(0, data['x'].shape[0]):


    rho[i] = 0. #initialize the density to 0

    '''
    there is no need to iterate a second time over the full array shape
    since the kernel is symmetric: ij = ji
    '''
    for j in range(0,i+1):

        dx = x[0,i]-x[0,j]
        dy = x[1,i]-x[1,j]
        dz = x[2,i]-x[2,j]
        dr = (dx**2+dy**2+dz**2)**0.5

        #calculate the density using the formula you saw in class
        rho[i] += m[j]*sigma/(pars[i_h]**3)*math.exp(-(dr)**2/(pars[i_h])**2)

        #to avoid double counting
        if i!=j:
                rho[j] += m[i]*sigma/(pars[i_h]**3)*math.exp(-(dr)**2/(pars[i_h])**2)


machs = []