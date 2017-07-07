import numpy as np
import math
import matplotlib.pyplot as plt


traject = np.loadtxt("T_2-ro_0.1.csv",delimiter=",")
n=4
N=n*n*n
nsamples=500


#diccionario de matrices Rij[time]

Rij={}

#recorro los tiempos
for time in range(0,nsamples):
 rij=np.zeros((N,N))    
#lleno la matriz rij del tiempo time
 for i in range(0,N):
  for j in range (0,i):
   dx=traject[i+N*time,1]-traject[j+N*time,1]
   dy=traject[i+N*time,2]-traject[j+N*time,2]
   dz=traject[i+N*time,3]-traject[j+N*time,3]
   rij[i,j]=np.sqrt((dx*dx)+(dy*dy)+(dz*dz)) #elemento ij de la matriz
 Rij[time]=rij #guardo la matriz


#matriz lindeman de pares:
 L={}
#En la posicion ij de L[time] va a contener [rij(time)-rij(0)/rij(0)]^2

#recorro los tiempos
for time in range(0,nsamples):
 lij=np.zeros((N,N))    
#lleno la matriz rij del tiempo time
 for i in range(0,N):
  for j in range (0,i):
   lij[i,j]=((Rij[time][i,j]-Rij[0][i,j])/(Rij[0][i,j]))**2#elemento ij de la matriz
 L[time]=lij #guardo la matriz

#Matriz de promedios
Prom_Rij={}

#lleno Prom_Rij con matrices nulas de NxN
for time in range(0,nsamples):
 promrij=np.zeros((N,N))    
 for i in range(0,N):
  for j in range (0,i):
   promrij[i,j]=0 #elemento ij de la matriz
 Prom_Rij[time]=promrij #guardo la matriz

#lleno Prom_Rij
suma=0
sumacuad=0
for i in range(0,N):
 for j in range(0,i):
  for time in range(0,nsamples):
    suma=suma+L[time][i,j]
    promedio=float(suma)/float(1+time)
    Prom_Rij[time][i,j]=promedio
  suma=0
  sumacuad=0 

#Ahora calculo el lindeman del tiempo time del solido
liendemang=[]

l=0
for time in range(0,nsamples): 
 for i in range(0,N):
  for j in range(0,i):
    l=l+Prom_Rij[time][i,j]                        #sumo
 liendemang.append([time,np.sqrt(l*2/(N*(N-1)))])  #promedio en el total de intercc N(N-1)/2 y raiz cuadrada
 l=0

#grafico
for time in range(0,nsamples):
 plt.plot(liendemang[time][0],liendemang[time][1],'ro')
plt.show()



