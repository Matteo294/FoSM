from matplotlib import pyplot as plt
from pandas import read_csv

data = read_csv("data.csv")
plt.plot(data['x1'], data['y1'], label='Planet', lw=2.0, color='blue')
plt.plot(data['x2'], data['y2'], label='Sun', lw=2.0, color='yellow')
ax=plt.gca()
ax.tick_params(direction='in')
ax.set_facecolor('black')
plt.xlabel('x [AU]', fontsize=14)
plt.ylabel('y [AU]', fontsize=14)
plt.legend(fontsize=14)
plt.show()

plt.plot(range(len(data['deltaE'])), data['deltaE'], color='black', lw=2.0)
plt.ylim([-0.001, 0.001])
plt.xlabel(r'Time   [$\times 2.98 \cdot 10^4 \, s$]', fontsize=14)
plt.ylabel(r'$\Delta E / E_0$', fontsize=14)
ax=plt.gca()
ax.tick_params(direction='in')
plt.show()