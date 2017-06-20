#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../../src/master.h"

double fuerza(double x);

int main()
{
  printf("--->Data\n");
  int numpoints = 10;
  double L = 100;
  double tabla[2*numpoints];
  printf("L: %f, Cant ptos tabla: %d\n",L,numpoints);
  printf("Total de particulas: 7\n");
  printf("Cercanas al borde de la caja: 6\n");
  printf("Centro de la caja: 1\n");
  make_table(&fuerza,numpoints,L,tabla);
  int i;
  printf("--->Tabla de fuerzas\n");
  printf("r\t\tf(r)\n");
  for(i=0; i<numpoints; i++)
  {
    printf("%f\t%.1f\n",tabla[i],tabla[i+numpoints]);
  }
  // Creo 7 particulas. 6 Muy cercanas al borde de la caja
  // y una en el centro.
  long int N = 7;
  struct part molec[N];
  double delta = 1.0;
  printf("--->Delta al techo o plano: %f\n",delta);
  printf("--->Indices\n");
  //Centro de la caja
  printf("CentroL/2 idx:0\n");
  molec[0].x = L/2;
  molec[0].y = L/2;
  molec[0].z = L/2;

  //Casi plano y-z
  printf("Plano y-z idx:1\n");
  molec[1].x = delta;
  molec[1].y = L/2;
  molec[1].z = L/2;

  //Casi plano x-y
  printf("Plano x-y idx:2\n");
  molec[2].x = L/2;
  molec[2].y = L/2;
  molec[2].z = delta;

  //Casi plano x-z
  printf("Plano x-z idx:3\n");
  molec[3].x = L/2;
  molec[3].y = delta;
  molec[3].z = L/2;

  //Casi techo x-y
  printf("Techo x-y idx:4\n");
  molec[4].x = L/2;
  molec[4].y = L/2;
  molec[4].z = L - delta;

  //Casi techo x-z
  printf("Techo x-z idx:5\n");
  molec[5].x = L/2;
  molec[5].y = L - delta;
  molec[5].z = L/2;

  //Casi techo z-y
  printf("Techo z-y idx:6\n");
  molec[6].x = L - delta;
  molec[6].y = L/2;
  molec[6].z = L/2;

  //Velocidades y fuerzas iniciales
  for(i=0;i<N;i++)
  {
    molec[i].vx = 0;
    molec[i].vy = 0;
    molec[i].vz = 0;

    molec[i].fx = 0;
    molec[i].fy = 0;
    molec[i].fz = 0;
  }

  // Evaluacion de fuerzas. La central debe tener fuerza nula y el
  // resto la equivalente a un r de 2.
  eval_f(molec,N,L,tabla,numpoints);
  printf("--->Fuerzas\n");
  printf("Particula\tfx\tfy\tfz\n");
  for(i=0;i<N;i++)
  {
    printf("%d\t\t%.2f\t%.2f\t%.2f\n",i,molec[i].fx,molec[i].fy,molec[i].fz);
  }

  return 0;
}

double fuerza(double x){
  if(x>20){
    return 0;
  }else{
    return 1;
  }
}
