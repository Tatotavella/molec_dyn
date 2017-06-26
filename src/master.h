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
	double ep; // Energia potencial
	double ec; // Energia cinetica
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
/*
 * @brief Escribe las fuerzas y potenciales en el struct @p *molec utilizando la
 * @brief discretizacion provista en @p *tabla de fuerza y potencial.
 * @brief Se utilizan condiciones periodicas de contorno.
 *
 * Se setean fuerzas dentro del struct @p *molec a 0.
 * La fuerza se asigna en la direccion que une a
 * las particulas.
 *
 * @date 19 Jun 2017
 * @author Franco Tavella
 */

int make_table(double(*funcion)(double), int numpoints, double L, double *tabla);
// Realiza una discretizacion de "func" en nptos con paso dr=L/nptos y los escribe
// en tabla

double funcion_LJ(double r);//Devuelve el potencial para un cierto valor de r. Realiza el corte en el potencial a partir de rc.

double funcion_fuerza(double r);//Devuelve la fuerza para un cierto valor de r. Realiza el corte en la fuerza a partir de rc.


int new_pos(struct part *past, struct part *future, long int N, double L, double h);
/*
 * @brief Escribe posiciones futuras en el @p *future a partir de aquellas del
 * @brief instante  @p h anterior @p *past. Aplica condiciones periodicas de
 * @brief contorno si la nueva posicion excede el tamaño @p L de la caja
 *
 * Para actualizar las posiciones se utiliza un desarrollo de Taylor a
 * segundo orden en el paso temporal @p h y las posiciones, velocidades
 * y fuerzas evaluadas en el instante pasado, @p *past.
 * Esta funcion se utiliza en una parte del algoritmo de Verlet en velocidades.
 *
 * @date 19 Jun 2017
 * @author Franco Tavella
 */

int new_vel(struct part *past, struct part *future, long int N, double L, double h);
/*
 * @brief Escribe velocidades futuras en el @p *future a partir de aquellas del
 * @brief instante  @p h anterior @p *past. Se utiliza la actualizacion
 * @brief correspondiente al algoritmo de Verlet en velocidades.
 *
 * Para actualizar las posiciones se utiliza un desarrollo de Taylor a
 * primer orden en medio paso temporal @p h y en un paso entero. De
 * esta manera la velocidad en un instante posterior se puede hallar
 * utilizando un promediado de las fuerzas en el instante pasado
 * y futuro. Antes de utilizar esta funcion se deben actualizar las
 * fuerzas en el struc @p *future.
 *
 * @date 19 Jun 2017
 * @author Franco Tavella
 */


double ord_verlet(struct part *molec, long int N, double L);
//Devuelve el parametro de ordenamiento geometrico de Verlet

int rescale_vel(struct part *molec, long int N, double T_old, double T_new);
//Escribe velocidades re-escaleadas en el srtuct


/* Faltan:
 * Funcion H de Boltzmann
 * Distribucion Radial
 * Re-escaleo de velocidades
 * y muchas mas
 */
