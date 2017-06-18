#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "master.h"

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

int init_rv(struct part *molec, long int N, double (*func)(double,double), double L, double T)
{
    int n = ceil(cbrt(N));
    double a = L/n;

    double vx_avg = 0;
    double vy_avg = 0;
    double vz_avg = 0;

    for (long int i = 0; i < N; i++) {
        long int l = i / (n * n);
        long int k = (i - l * n * n) / n;
        long int j = (i - l * n * n) % n;
        molec[i].x = j*a + a/2;
        molec[i].y = k*a + a/2;
        molec[i].z = l*a + a/2;

        molec[i].vx = (*func)(molec[i].m, T);
        molec[i].vy = (*func)(molec[i].m, T);
        molec[i].vz = (*func)(molec[i].m, T);
        vx_avg += molec[i].vx;
        vy_avg += molec[i].vy;
        vz_avg += molec[i].vz;
    }

    vx_avg = vx_avg / N;
    vy_avg = vy_avg / N;
    vz_avg = vz_avg / N;

    for (long int i = 0; i < N; i++) {
        molec[i].vx -= vx_avg;
        molec[i].vy -= vy_avg;
        molec[i].vz -= vz_avg;
    }

    return 0;
}
