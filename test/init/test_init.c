#include <stdlib.h>
#include <stdio.h>
#include "master.h"
#include "distribuciones.h"

int main()
{
    int n = 12;
    long int N = n*n*n;
    struct part *molec = malloc(N * sizeof(*molec));
    for (long int i = 0; i < N; i++) {
        molec[i].m = 1.0;
    }
    double L = 4.0;
    double T = 10.0;
    init_rv(molec, N, &maxwell_boltzmann, L, T);

    FILE *file = fopen("init_out.csv", "w");
    for (long int i = 0; i < N; i++) {
        fprintf(file, "%f,%f,%f,%f,%f,%f\n", molec[i].x, molec[i].y, molec[i].z, molec[i].vx, molec[i].vy, molec[i].vz);
    }
    fclose(file);
    free(molec);
    return 0;
}
