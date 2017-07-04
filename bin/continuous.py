import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    root_path = sys.argv[1]
else:
    root_path = './'
m = np.loadtxt(root_path + '/continuous_data.csv', delimiter=',')

n = 7
L = 9.5
N = n**3
V = L**3
rho = (n / L)**3

I = np.arange(m.shape[0]) + 1

K = m[:, 0]
V = m[:, 1]
E = m[:, 2]
Pex = m[:, 3]

T = 2 * K / 3 / N
P1 = (Pex + 1) * rho * T
P2 = Pex + rho * T
PI = rho * T

plt.figure()
plt.title('Energia')
plt.plot(I, E, '-o', color='C0', markersize=2.0, label='$E$')
plt.plot(I, K, '-o', color='C1', markersize=2.0, label='$K$')
plt.plot(I, V, '-o', color='C2', markersize=2.0, label='$V$')
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Presion')
plt.plot(I, Pex, '-o', color='C0', markersize=2.0, label='$P_{ex}$')
plt.plot(I, P1, '-o', color='C1', markersize=2.0, label='$P$ (gll)')
plt.plot(I, P2, '-o', color='C2', markersize=2.0, label='$P$ (mne)')
plt.plot(I, PI, '-o', color='C3', markersize=2.0, label='$P$ (idg)')
plt.grid()
plt.legend()
plt.show()
