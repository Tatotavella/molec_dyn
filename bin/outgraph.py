import numpy as np
import math
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
ekin = m[:,8]
epot = m[:,9]


plt.figure(1)
plt.plot(step,vx,'r')
plt.plot(step,vy,'g')
plt.plot(step,vz,'b')
plt.title("Velocidades")

plt.figure(2)
plt.plot(step,fx,'r')
plt.plot(step,fy,'g')
plt.plot(step,fz,'b')
plt.title("Fuerzas")

plt.figure(3)
plt.plot(step,ordver,'r')
plt.title("Ordenamiento Verlet")

plt.figure(4)
plt.plot(step,ekin,'r')
plt.plot(step,epot,'g')
plt.plot(step,ekin+epot,'b')
plt.title("Energia")

#H de Boltzman
hb = np.loadtxt("HB.csv",delimiter=",")
histvel = np.loadtxt("histvel.csv",delimiter=",")
stephb=hb[:,0]
hboltz = hb[:, 1]
histbin=histvel[:,0]
histcuenta=histvel[:,1]

#Distribucion de Boltzman:
ro=0.8442; #densidad
T=1;     #temp
f=[]
pi=math.pi;
vel=np.arange(0,4,0.1)

for i in range (0,len(vel)):
 f.append(4*ro*pi*math.pow(vel[i],2)*math.exp(-(math.pow(vel[i],2))/(2*T))/(math.pow(2*pi*T,1.5)))

plt.figure(5)
binsize=0.4
plt.bar(histbin,histcuenta,binsize)
plt.plot(vel,f,'r',label="MB T="+"%.2f" % float(T)+" "+r'$\rho=$'+"%.2f" %float(ro))
plt.title("Histograma velocidades")
plt.legend()
         
plt.figure(6)
plt.plot(stephb,hboltz,'g')
plt.title("H-Boltzman")
plt.show()

