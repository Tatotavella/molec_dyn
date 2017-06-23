#include <stdlib.h>
#include <stdio.h>
#include "master.h"
#include "distribuciones.h"

int main()
{
    long int n_p = 5;
    long int N = n_p*n_p*n_p;
    double T = 5;
    double L = 50;
    int Nr = L/0.001; //Con rcut asi hay puntos donde importa
    int Niter = 100000;
    double h = 0.01;

    int particul = 0;

    struct part *past = malloc(N * sizeof(*past));
    struct part *future = malloc(N * sizeof(*future));
    for (int i = 0; i < N; i++) {
        past[i].m = 1;
        future[i].m = 1;
    }
    init_rv(past, N, &maxwell_boltzmann, L, T);

    double *F = malloc(2*Nr * sizeof(*F));
    double *V = malloc(2*Nr * sizeof(*V));
    make_table(&funcion_fuerza, Nr, L, F);
    make_table(&funcion_LJ, Nr, L, V);
    eval_f(past, N, L, F, Nr);

    double **vel_avg = malloc((Niter+1) * sizeof(*vel_avg));
    double **vel_var = malloc((Niter+1) * sizeof(*vel_var));
    double **frz_avg = malloc((Niter+1) * sizeof(*frz_avg));
    double **frz_var = malloc((Niter+1) * sizeof(*frz_var));

    printf("running thermalization ...\n");
    FILE *file = fopen("out.csv", "w");
    double *lambda = malloc((Niter+1) * sizeof(*lambda));
    lambda[0] = ord_verlet(past, N, L);
    struct part *tmp = past;
    for (int i = 1; i < Niter+1; i++) {
        new_pos(past, future, N, L, h);
        eval_f(future, N, L, F, Nr);
        new_vel(past, future, N, L, h);
        tmp = past;
        past = future;
        future = tmp;

        lambda[i] = ord_verlet(past, N, L);

        vel_avg[i] = malloc(3 * sizeof(double));
        vel_var[i] = malloc(3 * sizeof(double));
        frz_avg[i] = malloc(3 * sizeof(double));
        frz_var[i] = malloc(3 * sizeof(double));
        promvar_v(past, N, vel_avg[i], vel_var[i]);
        promvar_f(past, N, frz_avg[i], frz_var[i]);

        printf("done with %d / %d\n", i, Niter);
        /*

        fprintf(file, "%d,%f,%f,%f,%f,%f,%f,%f\n", i, lambda[i],
                vel_avg[i][0], vel_avg[i][1], vel_avg[i][2],
                frz_avg[i][0], frz_avg[i][1], frz_avg[i][2]);
        */
        /*
        fprintf(file, "%d,%f,%f,%f,%f,%f,%f,%f\n", i, lambda[i],
                past[particul].vx, past[particul].vy, past[particul].vz,
                past[particul].fx, past[particul].fy, past[particul].fz);
        */
        //Energia
        double ecin = 0;
        double epot = 0;
        for (int i = 0; i < N; i++) {
            ecin += past[i].ec;
            epot += past[i].ep;
        }
        fprintf(file, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i, lambda[i],
        past[particul].vx, past[particul].vy, past[particul].vz,
        past[particul].fx, past[particul].fy, past[particul].fz,ecin,epot/2);
    }

    /*
    FILE *file = fopen("out.csv", "w");
    for (int i = 0; i < Niter+1; i++) {
        fprintf(file, "%d,%f,%f,%f,%f,%f,%f,%f\n", i, lambda[i],
                vel_avg[i][0], vel_avg[i][1], vel_avg[i][2],
                frz_avg[i][0], frz_avg[i][1], frz_avg[i][2]);
    }
    */
    fclose(file);

    free(past);
    free(future);
    free(F);
    free(V);
    free(lambda);
    free(vel_avg);
    free(vel_var);
    free(frz_avg);
    free(frz_var);

    return 0;
}
