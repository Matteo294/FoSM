from pandas import read_csv
from matplotlib import pyplot as plt 

data = read_csv("data.csv")
n = range(len(data['K']))
pos_before = read_csv("positions_before.csv")
pos_after = read_csv("positions_after.csv")

plt.plot(n, data['K'])
plt.xlabel("step")
plt.ylabel("K")
plt.show()

plt.plot(n, data['V'])
plt.xlabel("step")
plt.ylabel("V")
plt.show()

plt.plot(n, data['E'])
plt.xlabel("step")
plt.ylabel("E")
plt.show()

plt.plot(n, data['T'])
plt.plot([0, len(data['K'])], [83.0, 83.0])
plt.xlabel("step")
plt.ylabel("T")
plt.show()

plt.plot(n, data['Eerr'])
plt.xlabel("step")
plt.ylabel("deltaE")
plt.show()

ax = plt.axes(projection='3d')
for particle in pos_before.to_numpy():
    ax.scatter3D(particle[0], particle[1], particle[2])
plt.show()

ax = plt.axes(projection='3d')
for particle in pos_after.to_numpy():
    ax.scatter3D(particle[0], particle[1], particle[2], color='black')
plt.show()