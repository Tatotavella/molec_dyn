#ifndef DISTRIBUCIONES_H
#define DISTRIBUCIONES_H

/*! Devuelve un numero aleatorio uniformemente distribuido en el intervalo [0,1] */
double randu01();

/*! Devuelve un numero aleatorio distribuido segun la distribucion normal

    Para obtener muestras de la distribucion normal esta funcion implementa el
    algoritmo de Box-Muller.

    @param mean valor medio de la dsitribucion normal.
    @param std desviacion estandard de la distribucion normal.
*/
double randn(double mean, double std);

/*! Devuelve una velocidad distribuida segun la distribucion de Maxwell-Boltzmann

    La velocidad esta distribuida segun la distribucion de Maxwell-Boltzmann
    para una componente del vector velocidad de una particula de masa m en un
    gas ideal a temperatura T.

    @param m masa de las particulas.
    @param T temperatura.
*/
double maxwell_boltzmann(double m, double T);

#endif /* DISTRIBUCIONES_H */
