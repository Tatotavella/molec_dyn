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




int boltzmann(struct part *molec, long int N);




