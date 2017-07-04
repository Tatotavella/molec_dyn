import sys
import numpy as np
import matplotlib.pyplot as plt

if len(sys.argv) > 1:
    root_path = sys.argv[1]
else:
    root_path = './'
M = np.loadtxt(root_path + '/radialdist_data.csv', delimiter=',')

r = M[:, 0]
g_ini = M[:, 1]
g_fin = M[:, 2]

plt.figure(1)
plt.axhline(1, linestyle='--', color='gray')
plt.plot(r, g_ini, label='ini')
plt.plot(r, g_fin, label='fin')
plt.grid()
plt.legend()

plt.show()
