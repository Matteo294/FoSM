from os import read
from pandas import read_csv
from matplotlib import pyplot as plt 
import numpy as np

data_time = read_csv("times.csv")
data_errs = read_csv("errors.csv")

plt.scatter(range(1, len(data_time['midpoint'])+1), np.log(data_time['midpoint']), label='Midpoint Method', color='blue')
plt.scatter(range(1, len(data_time['mc'])+1), np.log(data_time['mc']), label='Monte Carlo', color='springgreen')
plt.legend(fontsize=12)
plt.xlabel('Dimensions', fontsize=14)
plt.ylabel('log(time) [ns]', fontsize=14)
plt.savefig("time.eps", dpi=300)
plt.show()

plt.scatter(range(1, len(data_errs['midpoint'])+1), data_errs['midpoint'], color='blue', label='Midpoint Method')
plt.scatter(range(1, len(data_errs['mc'])+1), data_errs['mc'], label='Monte Carlo', color='springgreen')
plt.legend(fontsize=12)
plt.xlabel('Dimensions', fontsize=14)
plt.ylabel('Error', fontsize=14)
plt.savefig("error.eps", dpi=300)
plt.show()