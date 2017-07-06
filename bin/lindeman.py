
import numpy as np
import math
import matplotlib.pyplot as plt


traject = np.loadtxt("trajectories_1.000000_data.csv",delimiter=",")
n=4
N=n*n*n
nsamples=1000

#diccionario de matrices Rij[time]

Rij={}

#recorro los tiempos
for time in range(0,nsamples):
 rij=np.zeros((N,N))    
#lleno la matriz rij del tiempo time
 for i in range(0,N):
  for j in range (0,i):
   dx=traject[i+(N-1)*time,1]-traject[j+(N-1)*time,1]
   dy=traject[i+(N-1)*time,2]-traject[j+(N-1)*time,2]
   dz=traject[i+(N-1)*time,3]-traject[j+(N-1)*time,3]
   rij[i,j]=np.sqrt(dx*dx+dy*dy+dz*dz) #elemento ij de la matriz
 Rij[time]=rij #guardo la matriz


#Matriz de promedios
Prom_Rij={}

#lleno Prom_Rij con matrices nulas de NxN
for time in range(0,nsamples):
 promrij=np.zeros((N,N))    
 for i in range(0,N):
  for j in range (0,i):
   promrij[i,j]=0 #elemento ij de la matriz
 Prom_Rij[time]=promrij #guardo la matriz

#Matriz de promedios cuad
Prom_Rijcuad={}

#lleno Prom_Rijcuad con matrices nulas de NxN
for time in range(0,nsamples):
 promrij=np.zeros((N,N))    
 for i in range(0,N):
  for j in range (0,i):
   promrij[i,j]=0 #elemento ij de la matriz
 Prom_Rijcuad[time]=promrij #guardo la matriz


#lleno Prom_Rij y Prom_Rijcuad
suma=0
sumacuad=0
for i in range(0,N):
 for j in range(0,i):
  for time in range(0,nsamples):
    suma=suma+Rij[time][i,j]
    sumacuad=sumacuad+(Rij[time][i,j]*Rij[time][i,j])
    promedio=float(suma)/float(1+time)
    promediocuad=float(sumacuad)/float(1+time)
    Prom_Rij[time][i,j]=promedio
    Prom_Rijcuad[time][i,j]=promediocuad
  suma=0
  sumacuad=0


#Ahora calculo el lindeman del tiempo time
liendeman=[]

l=0
for time in range(0,nsamples): 
 for i in range(0,N):
  for j in range(0,i):
    l=l+np.sqrt(np.abs(Prom_Rijcuad[time][i,j]-(Prom_Rij[time][i,j]*Prom_Rij[time][i,j])))/(Prom_Rij[time][i,j])   
 liendeman.append([time,l*2/(N*(N-1))])  
 l=0  

#grafico
for time in range(0,nsamples):
 plt.plot(liendeman[time][0],liendeman[time][1],'go')  
plt.show()
