import sys
import re
import glob
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    root_path = sys.argv[1]
else:
    root_path = './'

with open(root_path + '/parameters.txt', 'r') as f:
    for line in f:
        if line.startswith(r'N:'):
            m = re.search(r'N:(\d+)', line)
            N = int(m.group(1))
        elif line.startswith(r'L:'):
            m = re.search(r'L:(.+)', line)
            L = float(m.group(1))
        elif line.startswith(r'nsamples:'):
            m = re.search(r'nsamples:(.+)', line)
            nsamples = int(m.group(1))

n = int(np.sqrt(N))
V = L**3
rho = (n / L)**3

files_full_path = root_path + '/trajectories_*_data.csv'
files_list = glob.glob(files_full_path)

idx = -1
filename = files_list[idx]
print(filename)

filename = sys.argv[2]

times = np.arange(nsamples)
x = np.zeros(N)
y = np.zeros(N)
z = np.zeros(N)
xold = np.zeros(N)
yold = np.zeros(N)
zold = np.zeros(N)
x0 = np.zeros(N)
y0 = np.zeros(N)
z0 = np.zeros(N)
Rij = np.zeros((N, N, nsamples))
Rij2 = np.zeros((N, N, nsamples))
time = 0
with open(filename, 'r') as f:
    for line in f:
        val = line.split(',')
        pid = int(val[0])
        x[pid-1] = float(val[1])
        y[pid-1] = float(val[2])
        z[pid-1] = float(val[3])
        if np.abs(x[pid-1] - xold[pid-1]) >= L/2 and time > 0:
            x0[pid-1] += L * np.sign(xold[pid-1] - x[pid-1])
        if np.abs(y[pid-1] - yold[pid-1]) >= L/2 and time > 0:
            y0[pid-1] += L * np.sign(yold[pid-1] - y[pid-1])
        if np.abs(z[pid-1] - zold[pid-1]) >= L/2 and time > 0:
            z0[pid-1] += L * np.sign(zold[pid-1] - y[pid-1])
        xold[pid-1] = x[pid-1]
        yold[pid-1] = y[pid-1]
        zold[pid-1] = z[pid-1]
        x[pid-1] += x0[pid-1]
        y[pid-1] += y0[pid-1]
        z[pid-1] += z0[pid-1]
        if pid == N:
            dx = np.subtract.outer(x, x)
            dy = np.subtract.outer(y, y)
            dz = np.subtract.outer(z, z)
            dr2 = np.reshape(dx**2 + dy**2 + dz**2, (N, N, 1))
            Rij[:, :, time:] += np.sqrt(dr2)
            Rij2[:, :, time:] += dr2
            time += 1

for time in times:
    Rij[:, :, time] /= (time + 1)
    Rij2[:, :, time] /= (time + 1)

# calculo coeficiente lindemann del tiempo time
lindemann = np.zeros(nsamples)
L = (np.sqrt(Rij2 - Rij**2) / Rij) / (N*(N-1))
L[np.isnan(L)] = 0
lindemann = np.sum(L, axis=(0, 1))

# guardar data
np.save(filename + '_lindemann', np.array([times, lindemann]))

# grafico
plt.figure(1)
plt.plot(times, lindemann, '-')
plt.grid()
plt.show()
