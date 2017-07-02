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

    //Dist radial
    int bins = 500;
    double hist[bins];
    double n_hist[bins];
    float Ls = 1;
    for (int i=0; i<bins; i++){ 
        hist[i] = 0;
    }
    for (int i=0; i<bins; i++){
        n_hist[i] = i*Ls*L/bins;
    }
    //
    
    srand(time(NULL));

    struct part *past = malloc(N * sizeof(*past));
    struct part *future = malloc(N * sizeof(*future));
    struct thermo_data *data = malloc((Niter+1) * sizeof(*data));

    for (int i = 0; i < N; i++) {
        past[i].m = 1;
        future[i].m = 1;
    }
    init_rv(past, N, &maxwell_boltzmann, L, T);
    //dist_radial(past,N,L,bins,hist,n_hist,Ls);
    double *VF = malloc(3*Nr * sizeof(*VF));
    make_table(&funcion_LJ, &funcion_fuerza, Nr, L, VF);
    eval_f(past, N, L, VF, Nr, data);

    double **vel_avg = malloc((Niter+1) * sizeof(*vel_avg));
    double **vel_var = malloc((Niter+1) * sizeof(*vel_var));
    double **frz_avg = malloc((Niter+1) * sizeof(*frz_avg));
    double **frz_var = malloc((Niter+1) * sizeof(*frz_var));

    printf("running thermalization ...\n");
    char *filename = malloc((strlen(outdir) + 20) * sizeof(*filename));
    sprintf(filename, "%s/out.csv", outdir);
    FILE *file = fopen(filename, "w");
    FILE *file2 = fopen("gr.csv", "w");
    double *lambda = malloc((Niter+1) * sizeof(*lambda));
    lambda[0] = ord_verlet(past, N, L);
    FILE *fileHB = fopen("HB.csv", "w");//para funcion HBoltzman
    for (int i = 1; i < Niter+1; i++) {
        evolution_step(&past, &future, N, VF, Nr, L, h, (data+i));

        lambda[i] = ord_verlet(past, N, L);

        
        dist_radial(past,N,L,bins,hist,n_hist,Ls);
        
        
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

        //Energia
        double e[3];
        system_energy(past, N, e);
        */
        fprintf(file, "%d,%f,%f,%f,%f,%f,%f,%f,%f,%f\n", i, lambda[i],
        past[particul].vx, past[particul].vy, past[particul].vz,
        past[particul].fx, past[particul].fy, past[particul].fz, data[i].K, data[i].V);

	//Funcion H de Boltzman
	 //Histograma de las velocidades en cada paso
          float a=0;
          float b=4;
          int n=N;//numero de datos que se le pasan a histograma es igual al numero de molec, una vel por molecula
          int numbin=10;
          double bin=(b-a)/numbin;
	  //armo vector de velocidades
	  double *vel=malloc(N*sizeof(double));
	  for (int k=0;k<N;k++)
	    {
	      vel[k]=sqrt((past[k].vx*past[k].vx)+(past[k].vy*past[k].vy)+(past[k].vz*past[k].vz));
	    }
		 
          float *hist=malloc(2*numbin*sizeof(float));
          histograma(vel, hist, n, a, b, numbin);
	  
	  //Imprimo el ultimo cuando ya termalizo
	  if (i==Niter-1)
	    {
	     FILE *filehistvel = fopen("histvel.csv", "w");//para funcion HBoltzman
	     for(int j=0;j<numbin;j++)
	      {
		fprintf(filehistvel,"%f,%f\n",hist[numbin+j],hist[j]);
	      }
	     fclose(filehistvel);
	    }
	    
	  double HB=HBoltzman(past,hist,numbin,bin); //calcula la H Boltzman del paso i
          fprintf(fileHB, "%d,%f\n", i, HB);
	  free(hist);
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
    fclose(fileHB);
    free(filename);

    long int vmd_iter = 2000;
    struct thermo_data *data_vmd = malloc(vmd_iter * sizeof(*data_vmd));

    FILE *filemolec = fopen("molecevol.xyz", "w");
    for (int i = 0; i < vmd_iter; i++) {
        evolution_step(&past, &future, N, VF, Nr, L, h, (data_vmd+i));

        fprintf(filemolec, "%d\nmoleculas\n", N);
        for (int j = 0; j < N; j++) {
            fprintf(filemolec, "%d %f %f %f\n", j+1, past[j].x, past[j].y, past[j].z);
        }
    }
    fclose(filemolec);

    // Dist radial
    for (int i=0; i < bins; i++)
    {
        hist[i] = hist[i]/Niter;
    }
    double pi= 3.14159265358979323846;
    for (int i=0; i<bins; i++){
        hist[i] = (L*L*L)*hist[i]/(4*pi*n_hist[i]*n_hist[i]*(n_hist[1]-n_hist[0])*N*N/2);
    }
    for (int i=0; i<bins; i++){
        //printf("%f    %f\n",n_hist[i], hist[i]);
        fprintf(file2, "%f    %f\n",n_hist[i], hist[i]);
    }
    //

    free(past);
    free(future);
    free(data);
    free(data_vmd);
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
    int nsat, Lsat, drsat, hsat, Tsat, nitersat, osat;
    nsat = Lsat = drsat = hsat = Tsat = nitersat = osat = 0;
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
