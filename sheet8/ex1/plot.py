from pandas import read_csv
from matplotlib import pyplot as plt 

data = read_csv("data.csv")

plt.plot(data['x'], data['gauss'])
for iter in data.index.values[1:]:
    print(iter)
    plt.plot(data['x'], data[iter])
plt.show()