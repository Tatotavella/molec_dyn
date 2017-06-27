#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "master.h"
#include "distribuciones.h"

int parse_options(char ** argv, int argc, int *n, double *L, double *dr,
                  double *h, double *T, int *niter, char **outdir);

int main(int argc, char **argv)
{
    int n, N, Niter, Nr;
    double L, dr, h, T;
    char *outdir;

    int ret = parse_options(argv, argc, &n, &L, &dr, &h, &T, &Niter, &outdir);
    if (ret) return 1;
    N = n*n*n;
    Nr = L/dr;

    int particul = 0;

    srand(time(NULL));

    struct part *past = malloc(N * sizeof(*past));
    struct part *future = malloc(N * sizeof(*future));
    for (int i = 0; i < N; i++) {
        past[i].m = 1;
        future[i].m = 1;
    }
    init_rv(past, N, &maxwell_boltzmann, L, T);

    double *VF = malloc(3*Nr * sizeof(*VF));
    make_table(&funcion_LJ, &funcion_fuerza, Nr, L, VF);
    eval_f(past, N, L, VF, Nr);

    double **vel_avg = malloc((Niter+1) * sizeof(*vel_avg));
    double **vel_var = malloc((Niter+1) * sizeof(*vel_var));
    double **frz_avg = malloc((Niter+1) * sizeof(*frz_avg));
    double **frz_var = malloc((Niter+1) * sizeof(*frz_var));

    printf("running thermalization ...\n");
    char *filename = malloc((strlen(outdir) + 20) * sizeof(*filename));
    sprintf(filename, "%s/out.csv", outdir);
    FILE *file = fopen("out.csv", "w");
    double *lambda = malloc((Niter+1) * sizeof(*lambda));
    lambda[0] = ord_verlet(past, N, L);
    struct part *tmp = past;
    for (int i = 1; i < Niter+1; i++) {
        new_pos(past, future, N, L, h);
        eval_f(future, N, L, VF, Nr);
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
    free(filename);

    FILE *filemolec = fopen("molecevol.xyz", "w");
    for (int i = 0; i < 2000; i++) {
        new_pos(past, future, N, L, h);
        eval_f(future, N, L, VF, Nr);
        new_vel(past, future, N, L, h);
        tmp = past;
        past = future;
        future = tmp;

        fprintf(filemolec, "%d\nmoleculas\n", N);
        for (int j = 0; j < N; j++) {
            fprintf(filemolec, "%d %f %f %f\n", j+1, past[j].x, past[j].y, past[j].z);
        }
    }
    fclose(filemolec);

    free(past);
    free(future);
    free(VF);
    free(lambda);
    free(vel_avg);
    free(vel_var);
    free(frz_avg);
    free(frz_var);

    return 0;
}

int parse_options(char ** argv, int argc, int *n, double *L, double *dr,
                  double *h, double *T, int *niter, char **outdir)
{
    const char *usage = "usage: app [-n nparticles] [-L boxsize] [-dr spacing] [-h timestep] [-T temperature] [-niter niter] [-o outdir]\n"
                        "  -n nparticles: linear number of particles (total particles = n*n*n)\n"
                        "  -L boxsize: size of simulation box\n"
                        "  -dr spacing: spacing of the discretization of the force and potential\n"
                        "  -h timestep: discrete timestep of numerical integration\n"
                        "  -T temperature: initial temperature of the system\n"
                        "  -niter niter: number of iterations\n"
                        "  -o outdir: output directory where generated files will be stored\n";
    int nsat, Lsat, drsat, hsat, Tsat, nitersat, osat = 0;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-n") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *n = atoi(argv[i+1]);
            i++;
            nsat = 1;
        } else if (strcmp(argv[i], "-L") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *L = atof(argv[i+1]);
            i++;
            Lsat = 1;
        } else if (strcmp(argv[i], "-dr") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *dr = atof(argv[i+1]);
            i++;
            drsat = 1;
        } else if (strcmp(argv[i], "-h") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *h = atof(argv[i+1]);
            i++;
            hsat = 1;
        } else if (strcmp(argv[i], "-T") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *T = atof(argv[i+1]);
            i++;
            Tsat = 1;
        } else if (strcmp(argv[i], "-niter") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *niter = atoi(argv[i+1]);
            i++;
            nitersat = 1;
        } else if (strcmp(argv[i], "-o") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *outdir = argv[i+1];
            i++;
            osat = 1;
        }
    }

    if (!(nsat && Lsat && drsat && hsat && Tsat && nitersat && osat)) {
        printf("missing required arguments.\n");
        printf("%s", usage);
        return 1;
    }

    return 0;
}

