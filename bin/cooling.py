import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    root_path = sys.argv[1]
else:
    root_path = './'
m = np.loadtxt(root_path + '/cooling_data.csv', delimiter=',')

n = 7
L = 9.5
N = n**3
V = L**3
rho = (n / L)**3

T = m[:, 0]
K = m[:, 1]
V = m[:, 3]
E = m[:, 5]
Pex = m[:, 7]
Kstd = np.sqrt(m[:, 2])
Vstd = np.sqrt(m[:, 4])
Estd = np.sqrt(m[:, 6])
Pexstd = np.sqrt(m[:, 8])

P1 = (Pex + 1) * rho * T
P2 = Pex + rho * T
PI = rho * T

Tep = 2 * K / 3 / N

alpha = 1.96

plt.figure()
plt.plot(T, Tep, '-o')
plt.grid()
plt.show()

plt.figure()
plt.title('Energia')
plt.plot(T, E, '-o', color='C0', markersize=2.0, label='$E$')
plt.plot(T, K, '-o', color='C1', markersize=2.0, label='$K$')
plt.plot(T, V, '-o', color='C2', markersize=2.0, label='$V$')
plt.plot(T, K + V, '--')
plt.fill_between(T, E + alpha*Estd, E - alpha*Estd, color='C0', alpha=0.8)
plt.fill_between(T, K + alpha*Kstd, K - alpha*Kstd, color='C1', alpha=0.8)
plt.fill_between(T, V + alpha*Vstd, V - alpha*Vstd, color='C2', alpha=0.8)
plt.grid()
plt.legend()
plt.show()

plt.figure()
plt.title('Presion')
plt.plot(T, Pex, '-o', color='C0', markersize=2.0, label='$P_{ex}$')
plt.plot(T, P1, '-o', color='C1', markersize=2.0, label='$P$ (gll)')
plt.plot(T, P2, '-o', color='C2', markersize=2.0, label='$P$ (mne)')
plt.plot(T, PI, '-o', color='C3', markersize=2.0, label='$P$ (idg)')
plt.fill_between(T, Pex + alpha*Pexstd, Pex - alpha*Pexstd, color='C0', alpha=0.8)
plt.grid()
plt.legend()
plt.show()
