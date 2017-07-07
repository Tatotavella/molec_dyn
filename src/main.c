#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include "master.h"
#include "distribuciones.h"

int parse_options(char ** argv, int argc, int *n, double *L, double *dr,
                  double *h, double *Ti, double *Tf, int *Tsteps,
                  int *nsamples, int *nsep, int *ntherm, char **outdir);

int radial_distribution_bins(double L, int nbins, double *bins);

int radial_distribution_normalization(double L, long int N, int nsamples, int nbins, double *bins, double *freq);

int main(int argc, char **argv)
{
    /* Parse options and setup all program configuration options */
    int n, N, Nr, Tsteps, nsamples, nsep, ntherm;
    double L, dr, h, Ti, Tf, dT;
    char *outdir, *filename_continuous_data, *filename_cooling_data, *filename_radialdist_data, *filename_trajectories_data, *filename_parameters;
    FILE *file_continuous_data, *file_cooling_data, *file_radialdist_data, *file_trajectories_data, *file_parameters;

    int ret = parse_options(argv, argc, &n, &L, &dr, &h, &Ti, &Tf, &Tsteps, &nsamples, &nsep, &ntherm, &outdir);
    if (ret) return 1;
    N = n*n*n;
    Nr = L/dr;
    dT = -(Ti - Tf)/(Tsteps - 1);

    filename_parameters = malloc((strlen(outdir) + 20) * sizeof(*filename_parameters));
    sprintf(filename_parameters, "%s/parameters.txt", outdir);
    file_parameters = fopen(filename_parameters, "w");
    fprintf(file_parameters, "N:%d\nL:%f\ndr:%f\nh:%f\nTi:%f\nTf:%f\nTsteps:%d\nnsamples:%d\nnsep:%d\nntherm:%d\n",
            N, L, dr, h, Ti, Tf, Tsteps, nsamples, nsep, ntherm);
    fclose(file_parameters);
    free(filename_parameters);

    filename_continuous_data = malloc((strlen(outdir) + 20) * sizeof(*filename_continuous_data));
    sprintf(filename_continuous_data, "%s/continuous_data.csv", outdir);
    file_continuous_data = fopen(filename_continuous_data, "w");

    filename_cooling_data = malloc((strlen(outdir) + 20) * sizeof(*filename_cooling_data));
    sprintf(filename_cooling_data, "%s/cooling_data.csv", outdir);
    file_cooling_data = fopen(filename_cooling_data, "w");

    filename_radialdist_data = malloc((strlen(outdir) + 20) * sizeof(*filename_radialdist_data));
    sprintf(filename_radialdist_data, "%s/radialdist_data.csv", outdir);
    file_radialdist_data = fopen(filename_radialdist_data, "w");

    filename_trajectories_data = malloc((strlen(outdir) + 40) * sizeof(*filename_trajectories_data));

    srand(time(NULL));

    /* Initialize system */
    struct part *past = malloc(N * sizeof(*past));
    struct part *future = malloc(N * sizeof(*future));
    struct thermo_data observables;
    for (int i = 0; i < N; i++) {
        past[i].m = 1;
        future[i].m = 1;
    }
    init_rv(past, N, &maxwell_boltzmann, L, Ti, &observables);

    double *VF = malloc(3*Nr * sizeof(*VF));
    make_table(&funcion_LJ, &funcion_fuerza, Nr, L, VF);
    eval_f(past, N, L, VF, Nr, &observables);
    observables.E = observables.K + observables.V;

	int nbins = 500;
	double *rdist_bins = malloc(nbins * sizeof(*rdist_bins));
	double *rdist_freq_ini = malloc(nbins * sizeof(*rdist_freq_ini));
	double *rdist_freq_fin = malloc(nbins * sizeof(*rdist_freq_fin));
	radial_distribution_bins(L, nbins, rdist_bins);
	for (int i = 0; i < nbins; i++) {
		rdist_freq_ini[i] = 0;
		rdist_freq_fin[i] = 0;
	}

    /* Perform temperature sweep */
    printf("performing temperature sweep ...\n");
    struct thermo_data observables_avg = { 0, 0, 0, 0 };
    struct thermo_data observables_var = { 0, 0, 0, 0 };
    double Tnew, Tprev;
    for (int i = 0; i < Tsteps; i++) {
        /* New step temperature */
        Tnew = Ti + dT*i;

        /* Create trajectories file */
        sprintf(filename_trajectories_data, "%s/trajectories_%f_data.csv", outdir, Tnew);
        file_trajectories_data = fopen(filename_trajectories_data, "w");

        /* Thermalization */
        for (int j = 0; j < ntherm; j++) {
            if (j % 10 == 0) {
                Tprev = 2 * observables.K / (3 * N);
                rescale(past, N, Tprev, Tnew);
            }
            evolution_step(&past, &future, N, VF, Nr, L, h, &observables);
            fprintf(file_continuous_data, "%f,%f,%f,%f\n", observables.K, observables.V, observables.E, observables.Pex);
        }

        for (int j = 0; j < ntherm/10; j++) {
            evolution_step(&past, &future, N, VF, Nr, L, h, &observables);
            fprintf(file_continuous_data, "%f,%f,%f,%f\n", observables.K, observables.V, observables.E, observables.Pex);
        }

        /* Measurement */
        observables_avg = (struct thermo_data) { 0, 0, 0, 0 };
        observables_var = (struct thermo_data) { 0, 0, 0, 0 };

        for (int j = 0; j < nsamples; j++) {
            for (int k = 0; k < nsep; k++) {
                evolution_step(&past, &future, N, VF, Nr, L, h, &observables);
                fprintf(file_continuous_data, "%f,%f,%f,%f\n", observables.K, observables.V, observables.E, observables.Pex);
				if (i == 0) {
					dist_radial(past, N, L, nbins, rdist_freq_ini);
				} else if (i == Tsteps-1) {
					dist_radial(past, N, L, nbins, rdist_freq_fin);
				}
            }

            observables_avg.K += observables.K;
            observables_avg.V += observables.V;
            observables_avg.E += observables.E;
            observables_avg.Pex += observables.Pex;

            observables_var.K += observables.K * observables.K;
            observables_var.V += observables.V * observables.V;
            observables_var.E += observables.E * observables.E;
            observables_var.Pex += observables.Pex * observables.Pex;

            for (long int k = 0; k < N; k++) {
                fprintf(file_trajectories_data, "%ld,%f,%f,%f,%f,%f,%f,%f,%f\n", k+1,
                        past[k].x, past[k].y, past[k].z,
                        past[k].vx, past[k].vy, past[k].vz,
                        past[k].ec, past[k].ep);
            }
        }

        observables_avg.K /= nsamples;
        observables_avg.V /= nsamples;
        observables_avg.E /= nsamples;
        observables_avg.Pex /= nsamples;

        observables_var.K = observables_var.K / nsamples - observables_avg.K * observables_avg.K;
        observables_var.V = observables_var.V / nsamples - observables_avg.V * observables_avg.V;
        observables_var.E = observables_var.E / nsamples - observables_avg.E * observables_avg.E;
        observables_var.Pex = observables_var.Pex / nsamples - observables_avg.Pex * observables_avg.Pex;

        fprintf(file_cooling_data, "%f,%f,%f,%f,%f,%f,%f,%f,%f\n", Tnew,
                observables_avg.K, observables_var.K, observables_avg.V,
                observables_var.V, observables_avg.E, observables_var.E,
                observables_avg.Pex, observables_var.Pex);

        fclose(file_trajectories_data);

        printf("done with %d / %d\n", i+1, Tsteps);
    }

	radial_distribution_normalization(L, N, nsamples*nsep, nbins, rdist_bins, rdist_freq_ini);
	radial_distribution_normalization(L, N, nsamples*nsep, nbins, rdist_bins, rdist_freq_fin);
	for (int i = 0; i < nbins; i++) {
		fprintf(file_radialdist_data, "%f,%f,%f\n", rdist_bins[i], rdist_freq_ini[i], rdist_freq_fin[i]);
	}

    fclose(file_continuous_data);
    free(filename_continuous_data);

    fclose(file_cooling_data);
    free(filename_cooling_data);

	fclose(file_radialdist_data);
	free(filename_radialdist_data);

    free(past);
    free(future);
    free(VF);
	free(rdist_bins);
	free(rdist_freq_ini);
	free(rdist_freq_fin);

    return 0;
}

int radial_distribution_bins(double L, int nbins, double *bins)
{
    for (int i = 0; i < nbins; i++){ 
        bins[i] = i * L / 2 / nbins;
    }

	return 0;
}

int radial_distribution_normalization(double L, long int N, int nsamples, int nbins, double *bins, double *freq)
{
	double Vol = L * L * L;
	double dr = bins[1] - bins[0];

    for (int i = 0; i < nbins; i++) {
        freq[i] = Vol * freq[i] / (4 * M_PI * bins[i] * bins[i] * dr * N * N / 2) / nsamples;
    }

    return 0;
}

int rescale_expdecay(struct part *molec, long int N, double Told, double Tnew, double decay, double h)
{
    double alpha = sqrt(1 + (h/decay)*(Tnew - Told)/Told);

    for (long int i = 0; i < N; i++) {
        molec[i].vx *= alpha;
        molec[i].vy *= alpha;
        molec[i].vz *= alpha;
    }

    return 0;
}

int parse_options(char ** argv, int argc, int *n, double *L, double *dr,
                  double *h, double *Ti, double *Tf, int *Tsteps,
                  int *nsamples, int *nsep, int *ntherm, char **outdir)
{
    const char *usage = "usage: app [-n nparticles] [-L boxsize] [-dr spacing] [-h timestep]"
                        " [-Ti temperature] [-Tf temperature] [-Tsteps steps] [-nsamples samples]"
                        " [-nsep separation] [-ntherm thermalization steps] [-o outdir]\n"
                        "  -n nparticles: linear number of particles (total particles = n*n*n)\n"
                        "  -L boxsize: size of simulation box\n"
                        "  -dr spacing: spacing of the discretization of the force and potential\n"
                        "  -h timestep: discrete timestep of numerical integration\n"
                        "  -Ti temperature: initial temperature of the system\n"
                        "  -Tf temperature: final temperature of the system\n"
                        "  -Tsteps nsteps: number of steps between initial and final temperature\n"
                        "  -nsamples samples: number of samples to take at each temperature\n"
                        "  -nsep separation: number of points between consecutive samples\n"
                        "  -ntherm thermalization steps: number of steps for thermalization\n"
                        "  -o outdir: output directory where generated files will be stored\n";
    int nsat, Lsat, drsat, hsat, Tisat, Tfsat, Tstepssat, nsamplessat, nsepsat, nthermsat, osat;
    nsat = Lsat = drsat = hsat = Tisat = Tfsat = Tstepssat = nsamplessat = nsepsat = nthermsat = osat = 0;
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
        } else if (strcmp(argv[i], "-Ti") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *Ti = atof(argv[i+1]);
            i++;
            Tisat = 1;
        } else if (strcmp(argv[i], "-Tf") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *Tf = atof(argv[i+1]);
            i++;
            Tfsat = 1;
        } else if (strcmp(argv[i], "-Tsteps") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *Tsteps = atoi(argv[i+1]);
            i++;
            Tstepssat = 1;
        } else if (strcmp(argv[i], "-nsamples") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *nsamples = atoi(argv[i+1]);
            i++;
            nsamplessat = 1;
        } else if (strcmp(argv[i], "-nsep") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *nsep = atoi(argv[i+1]);
            i++;
            nsepsat = 1;
        } else if (strcmp(argv[i], "-ntherm") == 0) {
            if (argc <= i+1) {
                printf("failed option parsing.\n");
                printf("%s", usage);
                return 1;
            }
            *ntherm = atoi(argv[i+1]);
            i++;
            nthermsat = 1;
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

    if (!(nsat && Lsat && drsat && hsat && Tisat && Tfsat && Tstepssat && nsamplessat && nsepsat && nthermsat && osat)) {
        printf("missing required arguments.\n");
        printf("%s", usage);
        return 1;
    }

    return 0;
}

