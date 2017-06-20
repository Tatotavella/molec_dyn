#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "../../src/master.h"

double lineal(double x);

int main()
{
  int numpoints = 10;
  double L = 100;
  double tabla[2*numpoints];
  make_table(&lineal,numpoints,L,tabla);
  int i;
  printf("Tabla\n");
  printf("r\t\tf(r)\n");
  for(i=0; i<numpoints; i++)
  {
    printf("%f\t%f\n",tabla[i],tabla[i+numpoints]);
  }
  return 0;
  // Creo 7 particulas. 6 Muy cercanas al borde de la caja
  // y una en el centro.
  struct part molec[7];
}

double lineal(double x){
  return 2*x;
}
