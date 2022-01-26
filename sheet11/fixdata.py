import pandas as pd 
import numpy as np
from matplotlib import pyplot as plt
from os import listdir
from os.path import isfile, join

with open('Riemann.txt') as f:
    lines = f.read()
x = []
rho = []
p = []
u = []
for l in lines.splitlines():
    x.append(l[0])
    rho.append(l[1])
    p.append(l[2])
    u.append(l[3])

onlyfiles = [f for f in listdir("data") if isfile(join("data", f))]

for i in range(27):
    data = pd.read_csv("data/data_" + str(i) + ".csv")

    plt.plot(x, u, label="true")
    plt.plot(data['x'], data['u'])
    plt.legend()
    plt.show()

    plt.plot(x, rho, label="true")
    plt.plot(data['x'], data['rho'])
    plt.legend()
    plt.show()

    plt.plot(x, p, label="true")
    plt.plot(data['x'], data['p'])
    plt.legend()
    plt.show()
