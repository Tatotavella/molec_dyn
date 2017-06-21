#!/bin/env python3
import subprocess
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def gauss(x, mu, sigma):
    return np.exp(-((x-mu)**2)/2/sigma**2)/np.sqrt(2*np.pi*sigma**2)


def mb(v, m, T):
    return (4*np.pi*v**2)*np.exp(-(m*v**2)/(2*T))*(m/(2*np.pi*T))**(3/2)

subprocess.run('./test_init.e')

M = np.loadtxt('init_out.csv', delimiter=',')
x = M[:, 0]
y = M[:, 1]
z = M[:, 2]
vx = M[:, 3]
vy = M[:, 4]
vz = M[:, 5]
v = np.sqrt(vx**2 + vy**2 + vz**2)

T = 10
m = 1
var = T/m
std = np.sqrt(T/m)

print('{},{}'.format(x.min(), x.max()))
print('{},{}'.format(y.min(), y.max()))
print('{},{}'.format(z.min(), z.max()))

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, z)#markersize=4.0)
plt.show()

t = np.linspace(vx.min(), vx.max(), 500)
plt.figure()
plt.title('vx')
plt.hist(vx, bins=10, normed=True)
plt.plot(t, gauss(t, 0, std), '-')
plt.show()

t = np.linspace(vy.min(), vy.max(), 500)
plt.figure()
plt.title('vy')
plt.hist(vy, bins=10, normed=True)
plt.plot(t, gauss(t, 0, std), '-')
plt.show()

t = np.linspace(vz.min(), vz.max(), 500)
plt.figure()
plt.title('vz')
plt.hist(vz, bins=10, normed=True)
plt.plot(t, gauss(t, 0, std), '-')
plt.show()

t = np.linspace(v.min(), v.max(), 500)
plt.figure()
plt.title('v')
plt.hist(v, bins=10, normed=True)
plt.plot(t, mb(t, m, T), '-')
plt.show()
