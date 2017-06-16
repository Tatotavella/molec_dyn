#include <stdio.h>
#include <stdlib.h>
#include <time.h>

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

int init_rv(struct *part, long int N, double (*func)(double,double), double L, double T);
// Recibe una funcion y asigna al azar las velocidades de la distribucion "func" que recive (v,T). 
// Ordena las particulas equiespaciadas.

int promvar_v(struct *part, long int N, double *prom, double *var);
int promvar_f(struct *part, long int N, double *prom, double *var);
// Calculan promedio y varianza y lo escriben en prom y var. v y f respectivamente
// En prom y var escriben los 3 valores promedio. x -> 0 ; y -> 1 ; z -> 2

int eval_f(struct *part, long int N, double L, double *tabla);
//Interpola las fuerzas de la tabla y llena las propiedades del struct con eso.

int make_table(double (*func)(double), int nptos, double L, double *tabla);
// Realiza una discretizacion de "func" en nptos con paso dr=L/nptos y los escribe 
// en tabla

int new_pos(struct *past, struct *future, long int N, double L, double h);
//Escribe las nuevas posiciones en "future" basandose en "past" y taylor a segundo
//orden en la posicion

int new_vel(struct *past, struct *future, long int N, double L, double h);
//Escribe las nuevas velocidades en "future" basandose en "past" y taylor a primer
//orden en la velocidad con promedio de fuerzas.


double ord_verlet(struct *part, long int N, double L);
//Devuelve el parametro de ordenamiento geometrico de Verlet

/* Faltan:
 * Funcion H de Boltzmann
 * Distribucion Radial
 * Re-escaleo de velocidades
 * y muchas mas
 */








 







