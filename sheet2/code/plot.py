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

plt.plot(data['t'], data['deltaE'])
plt.show()