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

int make_table(double(*funcion_LJ)(double), double(*funcion_fuerza)(double), int numpoints, double L, double *tabla);
// Realiza una discretizacion de "funcion_LJ" y de "funcion_fuerza" en nptos con paso dr=L/nptos y los escribe
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

/*! Realiza un paso temporal de la evolucion del sistema

    @param past puntero al array de particulas en el estado inicial del paso.
        Al regresar la funcion, contiene el estado nuevo del sistema.
    @param future puntero al array de particulas que se utilizara para el calculo del
        neuvo estado. Al regresar la funcion, contiene el estado inicial del
        sistema.
    @param N numero de particulas en el sistema.
    @param VF array con la discretizacion de la fuerza y del potencial. El formato de
        la tabla esta dado por el output de la funcion @make_table.
    @param Nr numero de puntos en la discretizacion de la fuerza y potencial.
    @param L tamano de la caja de simulacion.
    @param h tamano del paso temporal de la integracion numerica.
*/
int evolution_step(struct part **past, struct part **future, long int N, double *VF, int Nr, double L, double h);

/*! Calcula la energia del sistema

    @param molec puntero al array de particulas del sistema.
    @param N numero de particulas en el sistema.
    @param e puntero a array donde se guarda el valor de la energia. El array
        debe tener tamano 3 y al regresar la funcion la entrada 0 contiene la
        energia total del sistema, la entrada 1 la energia cinetica y la
        entrada 2 la energia potencial.
*/
int system_energy(struct part* molec, long int N, double *e);

double ord_verlet(struct part *molec, long int N, double L);
//Devuelve el parametro de ordenamiento geometrico de Verlet

int dist_radial(struct part *molec, long int N, double L, int bins, double hist[], double n_hist[],float Ls);
// Genera la distribución radial.

int rescale(struct part *molec, long int N, double T_old, double T_new);
// Genera el reescaleo de velocidades.

/* Faltan:
 * Funcion H de Boltzmann
 * Re-escaleo de velocidades
 * y muchas mas
 */
double HBoltzman(struct part *molec, float *hist, long int numbin, double bin);
    /*Esta funcion va a recibir a past, a hist y al numero de bines y devolver el valor de HBoltzmann en ese "tiempo".
    Toma el histograma de las velocidades sobre las velocidades y a partir de este calcula H.
    */


int histograma(double *vel, float *hist, int n, float a, float b, int numbin);
/*Esta funcion realiza un histograma, con numcol la cantidad de bines, y extremos a y b, y cantidad de datos n=N
    //y[0]...y[m-1] cuentas
    //y[m]...y[m+m-1] marca de clase
    //n: numero de datos
    //[a,b] intervalo ext inf y sup
    //numcol: numero de columnas
    */
