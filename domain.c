#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "domain.h"
#include "constants.h"

/******************************************************************************
 * Read simulation domain dimensions, flag for the systems to use.
 *
 * Output:
 *  simul_domain: the domain information.
 ******************************************************************************/
void read_domain(int my_id, domain *simul_domain)
{
    FILE *fp;
    char *buff;
    int msg;
    fp = fopen("init.dat", "r");
    buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Domain information") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "xmin = %lf\n", &simul_domain->xmin);
    if (msg != 1) {
        printf("Failed to read xmin.\n");
        exit(1);
    }
    msg = fscanf(fp, "xmax = %lf\n", &simul_domain->xmax);
    if (msg != 1) {
        printf("Failed to read xmax.\n");
        exit(1);
    }
    msg = fscanf(fp, "ymin = %lf\n", &simul_domain->ymin);
    if (msg != 1) {
        printf("Failed to read ymin.\n");
        exit(1);
    }
    msg = fscanf(fp, "ymax = %lf\n", &simul_domain->ymax);
    if (msg != 1) {
        printf("Failed to read ymax.\n");
        exit(1);
    }
    msg = fscanf(fp, "zmin = %lf\n", &simul_domain->zmin);
    if (msg != 1) {
        printf("Failed to read zmin.\n");
        exit(1);
    }
    msg = fscanf(fp, "zmax = %lf\n", &simul_domain->zmax);
    if (msg != 1) {
        printf("Failed to read zmax.\n");
        exit(1);
    }
    if (my_id == 0) {
        printf("=========== Simulation Domain ===========\n");
        printf("xmin, xmax = %f %f\n", simul_domain->xmin, simul_domain->xmax);
        printf("ymin, ymax = %f %f\n", simul_domain->ymin, simul_domain->ymax);
        printf("zmin, zmax = %f %f\n", simul_domain->zmin, simul_domain->zmax);
        printf("=========================================\n");
    }

    /* /1* which system to use *1/ */
    /* while (fgets(buff, LEN_MAX, fp) != NULL) { */
    /*     //puts(buff); */
    /*     if (strstr(buff, "forcefree") != NULL) { */
    /*         break; */
    /*     } */
    /* } */
    /* msg = fscanf(fp, "Flag for the system to use:%d\n", &isystem); */
    /* if (msg != 1) { */
    /*     printf("Failed to read the flat for the system to use.\n"); */
    /*     exit(1); */
    /* } */
    /* fclose(fp); */
    /* free(buff); */

    /* /1* Time domain information *1/ */
    /* tracking_time_method(&simul_domain.tmax); */
    /* simul_domain.tmin = 0.0; */

    /* /1* Special information for each system *1/ */
    /* if (isystem == 1) { */
    /*     read_wlcs(); */
    /* } */
    /* else if (isystem == 2){ */
    /*     double v0; */
    /*     int tot_grid; */
    /*     dim_vfield(&v0); */
    /*     get_param_ff(); */
    /*     tot_grid = simul_grid.nx*simul_grid.ny*simul_grid.nz; */
    /*     vfd_b = (struct vfields*)malloc(sizeof(struct vfields)*tot_grid); */
    /*     read_vfields(0, v0, vfd_b); */
    /*     if (imultiple == 1) { */
    /*         vfd_a = (struct vfields*)malloc(sizeof(struct vfields)*tot_grid); */
    /*         read_vfields(1, v0, vfd_a); */
    /*     } */
    /* } */

    /* /1* Print the corresponding system information. *1/ */
    /* if (my_id == 0) { */
    /*     if (isystem == 0) { */
    /*         printf("Using system for test case...\n"); */
    /*     } */
    /*     else if (isystem == 1) { */
    /*         printf("Using wire-loop current systems...\n"); */
    /*     } else if (isystem==2){ */
    /*         printf("Using force-free magnetic field...\n"); */
    /*     } else { */
    /*         printf("ERROR: the system doesn't exit.\n" ); */
    /*         printf("Please choose a system from: 0 for test," */
    /*                 " 1 for wlcs, 2 for force-free"); */
    /*         exit(0); */
    /*     } */
    /* } */
}

/******************************************************************************
 * Read from init.dat to get the total particle tracking time and the tracking
 * method to use.
 *
 * Output:
 *  tott: the total particle tracking time.
 *  iadaptive: flag for tracking method to use.
 *             0 for fixed step. 1 for adaptive step.
 ******************************************************************************/
/* void tracking_time_method(double *tott) */
/* { */
/*     char buff[LEN_MAX]; */
/*     FILE *fp; */
/*     int msg; */
/*     fp = fopen("init.dat", "r"); */
/*     while (fgets(buff, LEN_MAX, fp) != NULL) { */
/*         //puts(buff); */
/*         if (strstr(buff, "Particle tracking") != NULL) { */
/*             break; */
/*         } */
/*     } */
/*     msg = fscanf(fp, "Total tracking time (s):%lf\n", tott); */
/*     if (msg != 1) { */
/*         printf("Failed to read the total tracking time.\n"); */
/*         exit(1); */
/*     } */
/*     msg = fscanf(fp, "Tracking method to use: %d\n", &iadaptive); */
/*     if (msg != 1) { */
/*         printf("Failed to read tracking method to use.\n"); */
/*         exit(1); */
/*     } */
/*     if (my_id == 0) { */
/*         printf("Total particle tracking time: %lf\n", *tott); */
/*         if (iadaptive == 0) { */
/*             printf("Fixed time step is used.\n"); */
/*         } */
/*         else if (iadaptive == 1) { */
/*             printf("Adaptive time step is used.\n"); */
/*         } */
/*         else { */
/*             printf("ERROR: wrong flag for time step."); */
/*             exit(0); */
/*         } */
/*     } */
/*     fclose(fp); */
/* } */
