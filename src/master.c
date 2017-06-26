/****************************************************************************
 *                                                                          *
 *                                                                          *
 ****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include "master.h"

int make_table(double(*funcion)(double), int numpoints, double L, double *tabla)
       /*
	*	Esta función recibe una funcion y completa el array *tabla el cual es
	*       de tamaño 2*numpoints. Además calcula la separación dr de los valores de
        *       la tabla mediante L y numpoints.
	*       El array *tabla se llenara de la siguiente manera:
        *         r_o, r_1,...,r_(numpoints-1), f(r_o), f(r_1), ...,f(r_numpoints-1)
	*/
 {
  double dr=L/numpoints;
  for (int i=0;i<numpoints;i++)
    {
      tabla[i]=dr*(i+1);                     //fila 1: x
      tabla[i+numpoints]=funcion(tabla[i]);  //fila 2: funcion(x)
    }
  return 0;
 }

double funcion_LJ(double r)
       /*
	*	Esta función devuelve el potencial de LJ a una distancia r entre particulas.
	*	Además utiliza una interpolacion con spline S(r)=ar^3+br^2+cr+d de orden 3
        *       entre los puntos rc=2.5sigma y rcorte=3sigma pidiendo:
        *       S(rc)=LJ(rc), S(rcorte)=0, S'(rc)=LJ'(rc) y S'(rcorte)=0.
        *       Esto permite que al hacer el corte del potencial la fuerza sea continua.
        *       El potencial V entonces queda como una funcion partida:
        *                 LJ(r)   si         r<rc
        *        V(r) =    S(r)   si     rc <r<rcorte
        *                   0     si  rcorte<r
	*/
 {
  double V;
  double sigma=1;
  double rc=2.5*sigma;
  double rcorte=3*sigma;
  //Coeficientes spline orden 3: ar^3+br^2+cr+d
  double a=-0.105072;
  double b=0.827847;
  double c=-2.130131;
  double d=1.776719;

  if (r<rc)
    {
     V=4*(pow(1/r,12)-pow(1/r,6));
    }
  else if (r<rcorte)
    {
     V=a*(r*r*r)+b*(r*r)+(c*r)+d;
    }
  else
    {
     V=0;
    }

 //posible forma de evitar ifs: int escalon=0.5*(1-(int)((x-rc)/fabs(x-rc)));  //escalon 1---0 centrada en x=rc

  return V;
 }

double funcion_fuerza(double r)
       /*
	*	Esta función devuelve la fuerza a una distancia r entre particulas.
	*	Además utiliza la derivada del spline S'(r)=ar^2+br+c (Ver "funcion_LJ")
        *       Esto permite que al hacer el corte del potencial la fuerza sea continua.
        *       La fuerza f entonces queda como una funcion partida:
        *                  -1*LJ'(r)   si         r<rc
        *        f(r) =    -1*S'(r)   si     rc <r<rcorte
        *                   0         si  rcorte<r
	*/
 {
  double f;
  double sigma=1;
  double rc=2.5*sigma;
  double rcorte=3*sigma;
  //coeficientes de S'(r) derivada del spline hecho en "funcion_LJ":
  double a=-0.105072;
  double b=0.827847;
  double c=-2.130131;

  if (r<rc)
    {
     f=(4*(12*pow(1/r,13)-6*pow(1/r,7)));
    }
  else if (r<rcorte)
    {
     f=-(3*a*(r*r)+2*b*(r)+(c));
    }
  else
    {
     f=0;
    }

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

int new_pos(struct part *past, struct part *future, long int N, double L, double h){
 long int i;
 double x_new, y_new, z_new;
 for(i=0; i<N; i++)
 {
   //Calculo de las nuevas posiciones
   x_new = past[i].x + past[i].vx * h + past[i].fx * ((h*h)/(2*past[i].m));
   y_new = past[i].y + past[i].vy * h + past[i].fy * ((h*h)/(2*past[i].m));
   z_new = past[i].z + past[i].vz * h + past[i].fz * ((h*h)/(2*past[i].m));

   //Asignacion mediante condiciones periodicas de contorno
   future[i].x = x_new - L * floor(x_new/L);
   future[i].y = y_new - L * floor(y_new/L);
   future[i].z = z_new - L * floor(z_new/L);
 }
 return 0;
}


int new_vel(struct part *past, struct part *future, long int N, double L, double h){
 long int i;
 double vx_new, vy_new, vz_new;
 double ecin;
 for(i=0; i<N; i++)
 {
   //Calculo de las nuevas velocidades
   vx_new = past[i].vx + (past[i].fx + future[i].fx) * (h/(2* past[i].m));
   vy_new = past[i].vy + (past[i].fy + future[i].fy) * (h/(2* past[i].m));
   vz_new = past[i].vz + (past[i].fz + future[i].fz) * (h/(2* past[i].m));

   //Asignacion
   future[i].vx = vx_new;
   future[i].vy = vy_new;
   future[i].vz = vz_new;

   //Energia cinetica
   ecin = (vx_new*vx_new + vy_new*vy_new + vz_new*vz_new)/(2*future[i].m);
   future[i].ec = ecin;

 }
 return 0;
}

int eval_f(struct part *molec, long int N, double L, double *tabla, int numpoints){
 long int i,j;
 double pre_force,potential;
 double x_dir,y_dir,z_dir;
 double dx,dy,dz;
 double x_box,y_box,z_box;
 double dr = L / numpoints;
 int index;
 double r_part;
 //Fuerzas a 0. Potencial a 0
 for (i = 0; i < N; i++)
 {
     molec[i].fx = 0;
     molec[i].fy = 0;
     molec[i].fz = 0;
     molec[i].ep = 0;
 }
 //Recorro todos los pares de particulas
 for(i=0; i<N; i++)
 {
   for(j=0; j<i; j++)
   {
    //Fuerza sobre las particulas i y j
    //Interaccion con particulas reales e imagen en una caja de tamaño L/2.
    //Posicion de las interactuantes. Correcion por condiciones periodicas.
    x_box = molec[j].x + L * (int)(2*(molec[i].x - molec[j].x)/L);
    y_box = molec[j].y + L * (int)(2*(molec[i].y - molec[j].y)/L);
    z_box = molec[j].z + L * (int)(2*(molec[i].z - molec[j].z)/L);

    dx = molec[i].x - x_box;
    dy = molec[i].y - y_box;
    dz = molec[i].z - z_box;

    r_part = sqrt(dx*dx + dy*dy + dz*dz);

    //Indice en la tabla para esa distancia
    index = ceil(r_part / dr) - 1;
    //Modulo de la fuerza y potencial para esa distancia. Ver funcion make_table
    pre_force = tabla[index + numpoints];
    //potential = tabla[index + 2*numpoints];

    /*
    //Potencial de la tabla
    molec[i].ep += potential;
    molec[j].ep += potential;
    */

    //Direccion de la fuerza
    x_dir = 1 * dx/r_part;
    y_dir = 1 * dy/r_part;
    z_dir = 1 * dz/r_part;

    //Asignacion de la fuerza para la particula i y j
    molec[i].fx += pre_force * x_dir;
    molec[i].fy += pre_force * y_dir;
    molec[i].fz += pre_force * z_dir;

    molec[j].fx += -1 * pre_force * x_dir;
    molec[j].fy += -1 * pre_force * y_dir;
    molec[j].fz += -1 * pre_force * z_dir;


    //Potencial calculado
    molec[i].ep += funcion_LJ(r_part);
    molec[j].ep += funcion_LJ(r_part);

    /*
    if(r_part<3)
    {
      printf("Interaccion part: %ld con %ld\n",i,j);
      printf("Distancia: %f, fuerza %f\n",r_part,pre_force);
      printf("Velocidades part: %ld vx: %f, vy: %f, vz: %f\n",i,
              molec[i].vx,molec[i].vy,molec[i].vz);
      printf("Velocidades part: %ld vx: %f, vy: %f, vz: %f\n",i,
              molec[i].vx,molec[i].vy,molec[i].vz);
      printf("Velocidades part: %ld vx: %f, vy: %f, vz: %f\n",j,
                      molec[j].vx,molec[j].vy,molec[j].vz);
    }
    */
   }
 }
 return 0;
}
