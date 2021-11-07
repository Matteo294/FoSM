from matplotlib import pyplot as plt
from pandas import read_csv

data = read_csv("data.csv")

plt.plot(data['x1'], data['y1'])
plt.plot(data['x2'], data['y2'])
plt.show()

plt.plot(range(len(data['energy'])), data['energy'])
plt.show()

plt.plot(range(len(data['deltaE'])), data['deltaE'])
plt.show()