import numpy as np
import matplotlib.pyplot as plt

times_tree = np.loadtxt('./times_tree.dat')

times_direct = np.loadtxt('./times_direct.dat')

N = np.array((5000,10000,20000,40000))

# lin reg:

k1, d1 = np.polyfit(np.log(N),np.log(times_tree[:,1]),1)
#k2, d2 = np.polyfit(np.log(N),np.log(times_direct[:,1]),1)
k2, d2 = np.polyfit(np.log(N),np.log(times_direct[:]),1)


fig = plt.figure()

plt.plot(N,times_tree[:,1],'bo',label="tree")
plt.plot(N,N**k1*np.exp(d1),'b-',label=f"{k1:.2f}*log(N)+{d1:.2f}")
#plt.plot(N,times_direct[:,1],'xr',label="direct")
plt.plot(N,times_direct[:],'xr',label="direct")
plt.plot(N,N**k2*np.exp(d2),'r-',label=f"{k2:.2f}*log(N)+{d2:.2f}")
plt.title(r'$\theta$ = 0.4')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('log(N)')
plt.ylabel('log(t)')
plt.legend(loc='upper left')

plt.show()

N = 1e10 

# t = N^k*exp(d)

t_tree10 = N**k1*np.exp(d1)

# y / d / h / m / s
y_tree = np.floor(t_tree10/31536000)
r = t_tree10 - y_tree*31536000
d_tree = np.floor(r/86400)
r = r - d_tree*86400
h_tree = np.floor(r/3600)
r = r - h_tree*3600
m_tree = np.floor(r/60)
r = r - m_tree*60
s_tree = r

t_direct10 = N**k2*np.exp(d2)
y_direc = np.floor(t_direct10/31536000)
r = t_direct10 - y_direc*31536000
d_direc = np.floor(r/86400)
r = r - d_direc*86400
h_direc = np.floor(r/3600)
r = r - h_direc*3600
m_direc = np.floor(r/60)
r = r - m_direc*60
s_direc = r

print(f'expected runtime for {N:.1e} particles with tree based calculation: {t_tree10}s ^= {y_tree}years {d_tree}days {h_tree}hours {m_tree}minutes {s_tree}seconds\n')
print(f'expected runtime for {N:.1e} particles with direct summation: {t_direct10}s ^= {y_direc}years {d_direc}days {h_direc}hours {m_direc}minutes {s_direc}seconds\n')




