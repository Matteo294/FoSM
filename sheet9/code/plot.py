from pandas import read_csv
from matplotlib import pyplot as plt 
import sys

mode = int(sys.argv[1]) # 0 for microcanonical 1 for canonical

data = read_csv("data.csv")
n = range(len(data['K']))
pos_before = read_csv("positions_before.csv")
pos_after = read_csv("positions_after.csv")

plt.plot(n, data['K'], label='Kinetic', linewidth=1.8, color='mediumspringgreen')
plt.plot(n, data['V'], label='Potential', linewidth=1.8, color='red')
plt.plot(n, data['E'], label='Total', linewidth=1.8, color='royalblue')
plt.title("Average energy per particle", fontsize=14)
plt.xlabel("step", fontsize=12)
plt.ylabel(r"E/$\epsilon$", fontsize=12)
plt.legend(fontsize=12)
ax = plt.gca()
ax.tick_params(direction='in')
plt.savefig("energies.eps", dpi=400)
plt.show()

plt.plot(n, data['T'], linewidth=1.8, color='black')
#plt.plot([0, len(data['K'])], [83.0, 83.0])
plt.xlabel("step", fontsize=12)
plt.ylabel("T", fontsize=12)
plt.title("Average temperature of the system", fontsize=14)
ax = plt.gca()
ax.tick_params(direction='in')
plt.savefig("temperature.eps", dpi=400)
plt.show()

plt.plot(n, data['Eerr'], linewidth=1.8, color='black')
plt.xlabel("step", fontsize=12)
plt.ylabel("deltaE", fontsize=12)
plt.title(r"Total energy deviation $\frac{E-E_0}{E_0}$", fontsize=14)
ax = plt.gca()
ax.tick_params(direction='in')
plt.savefig("deviation.eps", dpi=400)
plt.show()

ax = plt.axes(projection='3d')
for particle in pos_before.to_numpy():
    ax.scatter3D(particle[0], particle[1], particle[2], color='black')
ax = plt.gca()
ax.tick_params(direction='in')
plt.savefig("initial_conf.eps", dpi=400)
plt.show()

ax = plt.axes(projection='3d')
for particle in pos_after.to_numpy():
    ax.scatter3D(particle[0], particle[1], particle[2], color='black')
ax = plt.gca()
ax.tick_params(direction='in')
plt.savefig("final_conf.eps", dpi=400)
plt.show()