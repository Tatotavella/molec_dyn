import numpy as np
import matplotlib.pyplot as plt
m = np.loadtxt("out.csv",delimiter=",")
step = m[:,0]
ordver = m[:,1]
vx = m[:,2]
vy = m[:,3]
vz = m[:,4]
fx = m[:,5]
fy = m[:,6]
fz = m[:,7]
plt.plot(step,vx,'r')
plt.plot(step,vy,'g')
plt.plot(step,vz,'b')
plt.title("Velocidades")
plt.show()
plt.plot(step,fx,'r')
plt.plot(step,fy,'g')
plt.plot(step,fz,'b')
plt.title("Fuerzas")
plt.show()
plt.plot(step,ordver,'r')
plt.title("Ordenamiento Verlet")
plt.show()
