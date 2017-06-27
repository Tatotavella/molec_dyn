#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

struct part
{
	/* Propiedades de las particulas en un instante
	 *
	 *
	 *
	 */
  double x; //Posiciones
  double y;
	double z;
	double vx; //Velocidades
	double vy;
	double vz;
	double fx; //Fuerzas
	double fy;
	double fz;
	double e_p; // Energia potencial
	double e_c; // Energia cinetica
	double m; //Masa
};

int init_rv(struct part *molec, long int N, double (*func)(double,double), double L, double T);
// Recibe una funcion y asigna al azar las velocidades de la distribucion "func" que recive (v,T).
// Ordena las particulas equiespaciadas.

int promvar_v(struct part *molec, long int N, double *prom, double *var);
int promvar_f(struct part *molec, long int N, double *prom, double *var);
// Calculan promedio y varianza y lo escriben en prom y var. v y f respectivamente
// En prom y var escriben los 3 valores promedio. x -> 0 ; y -> 1 ; z -> 2

int eval_f(struct part *molec, long int N, double L, double *tabla, int numpoints);
//Interpola las fuerzas de la tabla y llena las propiedades del struct con eso.

int make_table(double(*funcion_LJ)(double), double(*funcion_fuerza)(double), int numpoints, double L, double *tabla);
// Realiza una discretizacion de "funcion_LJ" y de "funcion_fuerza" en nptos con paso dr=L/nptos y los escribe
// en tabla

double funcion_LJ(double r);//Devuelve el potencial para un cierto valor de r. Realiza el corte en el potencial a partir de rc.

double funcion_fuerza(double r);//Devuelve la fuerza para un cierto valor de r. Realiza el corte en la fuerza a partir de rc.


int new_pos(struct part *past, struct part *future, long int N, double L, double h);
//Escribe las nuevas posiciones en "future" basandose en "past" y taylor a segundo
//orden en la posicion

int new_vel(struct part *past, struct part *future, long int N, double L, double h);
//Escribe las nuevas velocidades en "future" basandose en "past" y taylor a primer
//orden en la velocidad con promedio de fuerzas.


double ord_verlet(struct part *molec, long int N, double L);
//Devuelve el parametro de ordenamiento geometrico de Verlet

int dist_radial(struct part *molec, long int N, double L, int bins, int hist[]);
// Genera la distribución radial.

int rescale(struct part *molec, long int N, double T_old, double T_new);
// Genera el reescaleo de velocidades.

/* Faltan:
 * Funcion H de Boltzmann
 * Re-escaleo de velocidades
 * y muchas mas
 */
