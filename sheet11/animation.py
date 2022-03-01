import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join
import matplotlib.animation as animation

def animate_rho(i):
    density.set_data(data[i]['x'], data[i]['rho'])
    return density,

def animate_u(i):
    speed.set_data(data[i]['x'], data[i]['u'])
    return speed,

def animate_p(i):
    pressure.set_data(data[i]['x'], data[i]['p'])
    return pressure,

writervideo = animation.FFMpegWriter(fps=15)


with open('Riemann.txt') as f:
    lines = f.read()
x = []
rho = []
p = []
u = []
for l in lines.splitlines():
    values = l.split()
    x.append(float(values[0]))
    rho.append(float(values[1]))
    p.append(float(values[2]))
    u.append(float(values[3]))

onlyfiles = [f for f in listdir("data") if isfile(join("data", f))]

data = []
for i in range(len(onlyfiles)):
    data.append(pd.read_csv("data/data_" + str(i) + ".csv"))

print("Rendering video 1/3")
fig = plt.figure()
density, = plt.plot(data[0]['x'], data[0]['rho'], label='FVM', color='mediumspringgreen')
plt.plot(x, rho, label="analytical", color='royalblue')
rhoAnimation = animation.FuncAnimation(fig, animate_rho, frames=len(onlyfiles), interval=5, blit=True, repeat=False)
plt.legend(fontsize=12)
plt.title("Density", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel(r"$\rho$", fontsize=12)
ax = plt.gca()
ax.tick_params(direction='in')
rhoAnimation.save('density.mp4', writer=writervideo)
plt.close()

print("Rendering video 2/3")
fig = plt.figure()
pressure, = plt.plot(data[0]['x'], data[0]['p'], label='FVM', color='mediumspringgreen')
plt.plot(x, p, label="analytical", color='royalblue')
pAnimation = animation.FuncAnimation(fig, animate_p, frames=len(onlyfiles), interval=5, blit=True, repeat=False)
plt.legend(fontsize=12)
plt.title("Pressure", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel(r"$p$", fontsize=12)
ax = plt.gca()
ax.tick_params(direction='in')
pAnimation.save('pressure.mp4', writer=writervideo)
plt.close()

print("Rendering video 3/3")
fig = plt.figure()
speed, = plt.plot(data[0]['x'], data[0]['u'], label='FVM', color='mediumspringgreen')
plt.plot(x, u, label="analytical", color='royalblue')
uAnimation = animation.FuncAnimation(fig, animate_u, frames=len(onlyfiles), interval=5, blit=True, repeat=False)
plt.legend()
plt.legend(fontsize=12)
plt.title("Speed", fontsize=14)
plt.xlabel("x", fontsize=12)
plt.ylabel(r"u", fontsize=12)
ax = plt.gca()
ax.tick_params(direction='in')
uAnimation.save('speed.mp4', writer=writervideo)
plt.close()

print("Done")
