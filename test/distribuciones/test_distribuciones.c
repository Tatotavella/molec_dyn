#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <assert.h>
#include "distribuciones.h"

int test_randn(double mean, double std);
int test_maxwell_boltzmann(double m, double T);

int main()
{
    srand(time(NULL));

    printf("tests randn...\n");
    test_randn(0, 1);
    printf("\n");
    test_randn(1, 10);
    printf("\n");
    test_randn(1, 0.1);
    printf("\n");

    printf("tests maxwell_boltzmann...\n");
    test_maxwell_boltzmann(1, 10);
    printf("\n");
    test_maxwell_boltzmann(1, 1);
    printf("\n");
    test_maxwell_boltzmann(1, 0.1);
    printf("\n");
}

int test_randn(double mean, double std)
{
    unsigned long int N = 10000000;
    double *v = malloc(N * sizeof(*v));
    for (unsigned long int i = 0; i < N; i++) {
        v[i] = randn(mean, std);
    }

    double avg = 0;
    for (unsigned long int i = 0; i < N; i++) {
        avg += v[i];
    }
    avg /= N;

    double m2 = 0;
    double m3 = 0;
    double m4 = 0;
    double m5 = 0;
    double dlt;
    for (unsigned long int i = 0; i < N; i++) {
        dlt = v[i] - avg;
        m2 += dlt*dlt;
        m3 += dlt*dlt*dlt;
        m4 += dlt*dlt*dlt*dlt;
    }
    m2 /= N;
    m3 /= N;
    m4 /= N;
    double skew = m3 / pow(m2, 3.0/2.0);
    double kurt = m4 / pow(m2, 2);

	free(v);

    printf("avg: %f == %f\n", avg, mean);
    if (mean != 0.0) {
	    assert(fabs(avg - mean)/fabs(mean) < 5E-3);
    } else {
	    assert(fabs(avg - mean) < 5E-3);
    }
    printf("var: %f == %f\n", m2, std*std);
	assert(fabs(m2 - std*std)/(std*std) < 5E-3);
    printf("skew: %f == %f\n", skew, 0.0);
	assert(fabs(skew) < 5E-3);
    printf("kurt: %f == %f\n", kurt, 3.0);
	assert(fabs(kurt - 3)/3 < 5E-3);

    return 0;
}

int test_maxwell_boltzmann(double m, double T)
{
    unsigned long int N = 10000000;
    double *v = malloc(N * sizeof(*v));
    for (unsigned long int i = 0; i < N; i++) {
        v[i] = maxwell_boltzmann(m, T);
    }

    double avg = 0;
    for (unsigned long int i = 0; i < N; i++) {
        avg += v[i];
    }
    avg /= N;

    double m2 = 0;
    double m3 = 0;
    double m4 = 0;
    double m5 = 0;
    double dlt;
    for (unsigned long int i = 0; i < N; i++) {
        dlt = v[i] - avg;
        m2 += dlt*dlt;
        m3 += dlt*dlt*dlt;
        m4 += dlt*dlt*dlt*dlt;
    }
    m2 /= N;
    m3 /= N;
    m4 /= N;
    double skew = m3 / pow(m2, 3.0/2.0);
    double kurt = m4 / pow(m2, 2);

	free(v);

    printf("avg: %f == %f\n", avg, 0.0);
	assert(fabs(avg) < 5E-3);
    printf("var: %f == %f\n", m2, T/m);
	assert(fabs(m2 - T/m)/(T/m) < 5E-3);
    printf("skew: %f == %f\n", skew, 0.0);
	assert(fabs(skew) < 5E-3);
    printf("kurt: %f == %f\n", kurt, 3.0);
	assert(fabs(kurt - 3)/3 < 5E-3);

    return 0;
}
