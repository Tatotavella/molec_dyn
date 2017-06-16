#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "estructuras.h"

int boltzmann(struct part *molec, long int N)
{
	/*
	*	

	*/
	int i = 6;
	molec[i].x = 10.0;
	molec[i].y = 10.0;
	molec[i].z = 10.0;
	molec[i].vx = 10.0;
	molec[i].vy = 10.0;
	molec[i].vz = 10.0;	
	

	return 0;
}

//OTRA FUNCION QUE GENERE LA TABLE 
// Vamos a usar dos struct para tener dos fotos del sistema en
// un instante anterior y posterio y la idea es actualizar
// los estados nuevos y viejos. Para el paso siguiente
// se convierte el "nuevo" en viejo intercambiando
// los punteros.

int lj_table(double *e_table, double eps, double sg, double cut, double dr)
{
	/*@brief devuelve la tabla llena de energias en e_table
	 * 
	 * eps y sg son las contantes del potencial	
	 * cut es el valor donde se corta la interaccion
	 * 
	 * OJO DEFINIR LA FORMA DE CORTE Y CORRECCION DEL POTENCIAL
	 */
	int i = 6;
	molec[i].x = 10.0;
	molec[i].y = 10.0;
	molec[i].z = 10.0;
	molec[i].vx = 10.0;
	molec[i].vy = 10.0;
	molec[i].vz = 10.0;	
	

	return 0;
}



int lj(struct part *molec, long int N, double *e_table, int p_i)
{
	/*@brief devuelve la energia de interaccion de Lennard Jones de p_i con 
	 * el reso de las particulas en el sistema.
	 * eps y sg son las contantes del potencial	
	 * cut es el valor donde se corta la interaccion
	 * OJO NECESITA LOOK UP TABLE ----> Le entra
	 * OJO ACA VAN CONDICIONES PERIODICAS DE CONTORNO
	 * OJO TAL VEZ PUEDE INTERPOLAR A LOS VALORES DE POTENCIAL
	 */
	int i = 6;
	molec[i].x = 10.0;
	molec[i].y = 10.0;
	molec[i].z = 10.0;
	molec[i].vx = 10.0;
	molec[i].vy = 10.0;
	molec[i].vz = 10.0;	
	

	return 0;
}

int main(int argc, char **argv) 
{

	long int N = 10; //Cantidad de particulas
	double L = 5.0; //Largo de la caja
	
	L = L + 1;
	
	//struct particula molec[N];
	struct part *molec_o = malloc( N * sizeof (struct part));
	struct part *molec_n = malloc( N * sizeof (struct part));

	int k = 6;
	

	printf("Antes\n");
	printf("Part : %d , x : %f, y : %f , z : %f , vx : %f , vy : %f , vz : %f\n",k,molec[k].x,molec[k].y,molec[k].z,molec[k].vx,molec[k].vy,molec[k].vz);
	//Inicializar los valores

	int i;

 	for(i = 0; i < N; i++)
	{
    	molec[i].x = 0.0;
		molec[i].y = 0.0;
		molec[i].z = 0.0;
		molec[i].vx = 0.0;
		molec[i].vy = 0.0;
		molec[i].vz = 0.0;
	}

    printf("Despues\n");
	printf("Part : %d , x : %f, y : %f , z : %f , vx : %f , vy : %f , vz : %f\n",k,molec[k].x,molec[k].y,molec[k].z,molec[k].vx,molec[k].vy,molec[k].vz);

	boltzmann(molec,N);

	printf("Despues Boltz\n");
	printf("Part : %d , x : %f, y : %f , z : %f , vx : %f , vy : %f , vz : %f\n",k,molec[k].x,molec[k].y,molec[k].z,molec[k].vx,molec[k].vy,molec[k].vz);

	return 0;
}
