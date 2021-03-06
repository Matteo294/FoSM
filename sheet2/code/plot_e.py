from matplotlib import pyplot as plt 
import pandas as pd 
import numpy as np

data = pd.read_csv("data.csv")
data_newIC = pd.read_csv("data_newIC.csv")

data['phi1'] = data['phi1'] * (180/np.pi)
data['phi2'] = data['phi2'] * (180/np.pi)
data_newIC['phi1'] = data_newIC['phi1'] * (180/np.pi)
data_newIC['phi2'] = data_newIC['phi2'] * (180/np.pi)

plt.plot(data['phi1'][0], data['phi2'][0], '.', markersize=16, color='red')
plt.plot(data_newIC['phi1'][0], data_newIC['phi2'][0], '.', markersize=16, color='blue')
plt.plot(data['phi1'], data['phi2'],'k-', label=r'IC: $\phi_1 = 50°$ & $\phi_2 = -120°$')
plt.plot(data_newIC['phi1'], data_newIC['phi2'],'y-', label=r'IC: $\phi_1 = 51°$ & $\phi_2 = -120°$')
plt.legend()
plt.xlabel(r'$\phi_1$ (°)')
plt.ylabel(r'$\phi_2$ (°)')
plt.show()
plt.plot(data['phi1dot'][0], data['phi2dot'][0], '.', markersize=16, color='red')
plt.plot(data_newIC['phi1dot'][0], data_newIC['phi2dot'][0], '.', markersize=16, color='blue')
plt.plot(data['phi1dot'], data['phi2dot'],'k-', label=r'IC: $\phi_1 = 50°$ & $\phi_2 = -120°$')
plt.plot(data_newIC['phi1dot'], data_newIC['phi2dot'],'y-', label=r'IC: $\phi_1 = 51°$ & $\phi_2 = -120°$')
plt.legend()
plt.xlabel(r'$\dot{\phi}_1$ (°s⁻¹)')
plt.ylabel(r'$\dot{\phi}_2$ (°s⁻¹)')
plt.show()

#plt.plot(data['t'], data['deltaE'], color='red', linewidth=1.8)
#plt.xlabel('t', fontsize=16)
#plt.ylabel(r'$\epsilon(t)$', fontsize=16)
#plt.grid()
#plt.show()
