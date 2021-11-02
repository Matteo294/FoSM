from matplotlib import pyplot as plt 
import pandas as pd 

data = pd.read_csv("data.csv")

plt.plot(data['t'], data['phi1'], label='phi1')
plt.plot(data['t'], data['phi2'], label='phi2')
plt.show()
plt.plot(data['t'], data['phi1dot'], label='phi1dot')
plt.plot(data['t'], data['phi2dot'], label='phi2dot')
plt.legend()
plt.xlabel('t')
plt.show()

plt.plot(data['t'], data['deltaE'], color='red', linewidth=1.8)
plt.xlabel('t', fontsize=16)
plt.ylabel(r'$\epsilon(t)$', fontsize=16)
plt.grid()
plt.show()