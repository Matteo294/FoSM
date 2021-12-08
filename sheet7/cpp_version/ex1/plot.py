from pandas import read_csv
from matplotlib import pyplot as plt

data1 = read_csv("gen1_data.csv")
data2 = read_csv("gen2_data.csv")
data_RANDU = read_csv("RANDU_data.csv")

plt.scatter(data1['x'], data1['y'], color='blue')
plt.xlim([0.2, 0.201])
plt.ylim([0.3, 0.301])
plt.savefig("mt19937.eps", dpi=300)
plt.show()
plt.scatter(data2['x'], data2['y'], color='mediumspringgreen')
plt.xlim([0.2, 0.201])
plt.ylim([0.3, 0.301])
plt.savefig("builtin.eps", dpi=300)
plt.show()
plt.scatter(data_RANDU['x'], data_RANDU['y'], color='red')
plt.xlim([0.2, 0.201])
plt.ylim([0.3, 0.301])
plt.savefig("randu.eps", dpi=300)
plt.show()