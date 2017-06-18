#include <stdlib.h>
#include <math.h>

double randu01()
{
    return ((double)rand())/RAND_MAX;
}

double randn(double mean, double std)
{
    double U = randu01();
    double V = randu01();
    double r = sqrt(- 2 * log(U)) * cos(2 * M_PI * V);
    return std*r + mean;
}

double maxwell_boltzmann(double m, double T)
{
    double std = sqrt(T/m);
    return randn(0, std);
}
