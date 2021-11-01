from matplotlib import pyplot as plt 
import pandas as pd
from matplotlib.animation import FuncAnimation

global data 
data = pd.read_csv("animation.csv")

Ltot = 5
xmin = -Ltot-0.5
xmax = Ltot+0.5
ymin = -Ltot-0.5
ymax = Ltot+0.5

fig = plt.figure()
ax = plt.axes()
axes = fig.gca()
axes.axis("equal")
axes.set_xlim(xmin, xmax)
axes.set_ylim(ymin, ymax)
axes.set_axis_off()
pendulum1, = ax.plot([0, 2], [0, 2], color='black', linewidth=3)
pendulum2, = ax.plot([0, 2], [0, 2], color='black', linewidth=3)
balls, = ax.plot([0, 0], [0, 0], 'o', markersize=15, color='black')

def init():
    pendulum1.set_data([], [])
    pendulum2.set_data([], [])
    balls.set_data([], [])
    return pendulum1, pendulum2, balls
def animate(i):
    pendulum1.set_data([0, data['x1'].to_numpy()[i]], [0, data['y1'].to_numpy()[i]])
    pendulum2.set_data([data['x1'].to_numpy()[i], data['x2'].to_numpy()[i]], [data['y1'].to_numpy()[i], data['y2'].to_numpy()[i]])
    balls.set_data([data['x1'].to_numpy()[i], data['x2'].to_numpy()[i]], [data['y1'].to_numpy()[i], data['y2'].to_numpy()[i]])
    return pendulum1, pendulum2, balls

anim = FuncAnimation(fig, animate, init_func=init, blit=True, interval=20)
plt.show()

data = pd.read_csv("animation.csv")

'''
plt.plot(data['x1'], data['y1'])
plt.plot(data['x2'], data['y2'])
plt.plot([0, data['x1'].to_numpy()[-1]], [2, data['y1'].to_numpy()[-1]], '-', color='black', linewidth=1.8)
plt.plot([data['x1'].to_numpy()[-1], data['x2'].to_numpy()[-1]], [data['y1'].to_numpy()[-1], data['y2'].to_numpy()[-1]], '-', color='black', linewidth=1.8)
plt.ylim((0, 2))
plt.show()
'''