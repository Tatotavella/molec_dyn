#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "estructuras.h"
#include <math.h>

int make_table(double(*funcion)(double), int numpoints, double L, double *tabla)
 {
  double dr=L/numpoints;
  for (int i=0;i<numpoints;i++)
    {
      tabla[i]=dr*(i+1);                     //fila 1: x ---->dr*(i+1)para que no evalue en 0
      tabla[i+numpoints]=funcion(tabla[i]);  //fila 2: funcion(x)
    }
  return 0;
 }

double funcion_LJ(double r)
 {
  double V;
  double rc=2.4999;                              //de esta forma no tira error en fabs(x-rc) si llega a rc=2.5 (ver si se puede mejorar el escalon)
  int escalon=0.5*(1-(int)((r-rc)/fabs(r-rc)));  //escalon 1---0 centrada en x=rc
  V=(4*(pow(1/r,12)-pow(1/r,6))-4*(pow(1/rc,12)-pow(1/rc,6)))*(escalon); //potencial por la funcion escalon

  return V;
 }

double funcion_fuerza(double r)
 {
  double f;
  double rc=2.4999;                       //de esta forma no tira error en fabs(x-rc) si llega a rc=2.5 (ver si se puede mejorar el escalon)
  int escalon=0.5*(1-(int)((r-rc)/fabs(r-rc)));  //escalon 1---0 centrada en x=rc
  
  f=(4*(12*pow(1/r,13)-6*pow(1/r,7)))*(escalon); //derivada del potencial por la funcion escalon

  return f;
 }
