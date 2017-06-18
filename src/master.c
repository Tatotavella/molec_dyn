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


int promvar_v(struct part *molec, long int N, double *prom, double *var){

	/*
	*	Esta función calcula el promedio y la varianza para un struct de N particulas
	*	para las velocidades vx, vy y vz independientemente.
	*	El resultado final es asignado a los punteros *prom y *var respectivamente.
	*   Ejemplo:
	*		double prom_vec[3];
	*		double *prom;
	*		prom = &prom_vec[0];
	*/

	int i;
	// Promedio.
	*(prom+0) = 0, *(prom+1) = 0, *(prom+2) = 0; // Se inicializan los valores en 0.
	for (i=0; i<N; i++)
	{
		*(prom+0) = *(prom+0) + molec[i].vx;
		*(prom+1) = *(prom+1) + molec[i].vy;
		*(prom+2) = *(prom+2) + molec[i].vz;
	}	
	*(prom+0) = *(prom+0)/N, *(prom+1) = *(prom+1)/N, *(prom+2) = *(prom+2)/N;
	
	// Varianza.
	*(var+0) = 0, *(var+1) = 0, *(var+2) = 0; // Se inicializan los valores en 0.
	for (i=0; i<N; i++)
	{
		*(var+0) = *(var+0) + (molec[i].vx-*(prom+0))*(molec[i].vx-*(prom+0));
		*(var+1) = *(var+1) + (molec[i].vy-*(prom+1))*(molec[i].vy-*(prom+1));
		*(var+2) = *(var+2) + (molec[i].vz-*(prom+2))*(molec[i].vz-*(prom+2));
	}	
	*(var+0) = *(var+0)/(N-1), *(var+1) = *(var+1)/(N-1), *(var+2) = *(var+2)/(N-1);

	return 0;

} 


int promvar_f(struct part *molec, long int N, double *prom, double *var){

	/*
	*	Esta función calcula el promedio y la varianza para un struct de N particulas
	*	para las fuerzas fx, fy y fz independientemente.
	*	El resultado final es asignado a los punteros *prom y *var respectivamente.
	*	Ejemplo:
	*		double prom_vec[3];
	*		double *prom;
	*		prom = &prom_vec[0];
	*/

	int i;
	// Promedio.
	*(prom+0) = 0, *(prom+1) = 0, *(prom+2) = 0; // Se inicializan los valores en 0.
	for (i=0; i<N; i++)
	{
		*(prom+0) = *(prom+0) + molec[i].fx;
		*(prom+1) = *(prom+1) + molec[i].fy;
		*(prom+2) = *(prom+2) + molec[i].fz;
	}	
	*(prom+0) = *(prom+0)/N, *(prom+1) = *(prom+1)/N, *(prom+2) = *(prom+2)/N;
	
	// Varianza.
	*(var+0) = 0, *(var+1) = 0, *(var+2) = 0; // Se inicializan los valores en 0.
	for (i=0; i<N; i++)
	{
		*(var+0) = *(var+0) + (molec[i].fx-*(prom+0))*(molec[i].fx-*(prom+0));
		*(var+1) = *(var+1) + (molec[i].fy-*(prom+1))*(molec[i].fy-*(prom+1));
		*(var+2) = *(var+2) + (molec[i].fz-*(prom+2))*(molec[i].fz-*(prom+2));
	}	
	*(var+0) = *(var+0)/(N-1), *(var+1) = *(var+1)/(N-1), *(var+2) = *(var+2)/(N-1);

	return 0;

}


double ord_verlet(struct part *molec, long int N, double L){

	/*
	*	Esta función realiza el calculo del ordenamiento de Verlet.
	*	Esta función devuelve 1 cuando el sistema se encuentra puramente
	* 	ordenado y 0 en el caso contrario.
	*   Importante: Para su correcta utilización y por la forma en que se
	*   define a (distancia entre moleculas), es importante que N sea un numero
	*   tal que su raiz cubica sea un número entero.
	*/

	double lamda_x=0,lamda_y=0,lamda_z=0,lamda;
	double pi = 3.14159265358979323846;
	double a  = L/pow(N,1.0/3.0);

	int i;
	for (i=0; i<N; i++)
	{
		lamda_x = lamda_x + cos((2*pi/a)*(molec[i].x-a/2));
		lamda_y = lamda_y + cos((2*pi/a)*(molec[i].y-a/2));
		lamda_z = lamda_z + cos((2*pi/a)*(molec[i].z-a/2));
	}
	lamda_x=lamda_x/N, lamda_y=lamda_y/N, lamda_z=lamda_z/N;
	lamda = (lamda_x+lamda_y+lamda_z)/3;

	return lamda;

}
