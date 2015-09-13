#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "domain.h"
#include "constants.h"
#include "wlcs.h"

void get_spatial_domain(int mpi_rank, FILE *fp, domain *simul_domain);
int get_system_type(int mpi_rank, FILE *fp);
int get_tracking_time_method(int mpi_rank, FILE *fp, double *tott,
        int *tracking_method);
int get_system_info(int mpi_rank, int system_type, char *config_file_name);

/******************************************************************************
 * Read simulation domain dimensions, flag for the systems to use.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  config_file_name: the configuration file name.
 *
 * Output:
 *  simul_domain: the domain information.
 ******************************************************************************/
int read_domain(int mpi_rank, char *config_file_name, domain *simul_domain)
{
    FILE *fp;
    int system_type, err, tracking_method;
    double tott;
    fp = fopen(config_file_name, "r");
    get_spatial_domain(mpi_rank, fp, simul_domain);
    system_type = get_system_type(mpi_rank, fp);
    if (system_type < 0) {
        fclose(fp);
        return -1;
    }
    err = get_tracking_time_method(mpi_rank, fp, &tott, &tracking_method);
    if (err < 0) {
        fclose(fp);
        return -1;
    }
    simul_domain->tmax = tott;
    simul_domain->tmin = 0.0;

    fclose(fp);

    err = get_system_info(mpi_rank, system_type, config_file_name);
    if (err < 0) {
        return -1;
    }

    return 0;
}

/******************************************************************************
 * Read simulation spatial domain information.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  fp: file handler.
 *
 * Output:
 *  simul_domain: the domain information.
 ******************************************************************************/
void get_spatial_domain(int mpi_rank, FILE *fp, domain *simul_domain)
{
    char *buff;
    int msg;
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
    if (mpi_rank == 0) {
        printf("=========== Simulation Domain ===========\n");
        printf("xmin, xmax = %f %f\n", simul_domain->xmin, simul_domain->xmax);
        printf("ymin, ymax = %f %f\n", simul_domain->ymin, simul_domain->ymax);
        printf("zmin, zmax = %f %f\n", simul_domain->zmin, simul_domain->zmax);
        printf("=========================================\n");
    }

    free(buff);
}


/******************************************************************************
 * Read from init.dat to get the total particle tracking time and the tracking
 * method to use.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  fp: file handler.
 *
 * Output:
 *  tott: the total particle tracking time.
 *  tracking_method: flag for tracking method to use.
 *             0 for fixed step. 1 for adaptive step.
 ******************************************************************************/
int get_tracking_time_method(int mpi_rank, FILE *fp, double *tott,
        int *tracking_method)
{
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    int msg;
    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Particle tracking") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "Total tracking time (s):%lf\n", tott);
    if (msg != 1) {
        printf("Failed to read the total tracking time.\n");
        exit(1);
    }
    fgets(buff, LEN_MAX, fp); // Skip the description
    msg = fscanf(fp, "Tracking method to use: %d\n", tracking_method);
    if (msg != 1) {
        printf("Failed to read tracking method to use.\n");
        exit(1);
    }
    if (mpi_rank == 0) {
        printf("Total particle tracking time: %lf\n", *tott);
    }
    switch (*tracking_method) {
        case 0:
            if (mpi_rank == 0) {
                printf("Fixed time step is used.\n");
            }
            break;
        case 1:
            if (mpi_rank == 0) {
                printf("Adaptive time step is used.\n");
            }
            break;
        default:
            if (mpi_rank == 0) {
                printf("ERROR: wrong flag for time step.");
            }
            return -1;
    }
    free(buff);
    return 0;
}

/******************************************************************************
 * Read simulation domain dimensions, flag for the systems to use.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  fp: file handler.
 *
 * Return:
 *  system_type: the system type.
 ******************************************************************************/
int get_system_type(int mpi_rank, FILE *fp)
{
    int system_type, msg;
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    /* which system to use */
    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "forcefree") != NULL) {
            break;
        }
    }
    system_type = 0;
    msg = fscanf(fp, "Flag for the system to use:%d\n", &system_type);
    if (msg != 1) {
        printf("Failed to read the flat for the system to use.\n");
        exit(1);
    }
    switch(system_type) {
        case 0:
            if (mpi_rank == 0) {
                printf("The system is a test system.\n");
            }
            break;
        case 1:
            if (mpi_rank == 0) {
                printf("The system is a wire-loop current system.\n");
            }
            break;
        case 2:
            if (mpi_rank == 0) {
                printf("The system is a force-free magnetic field system.\n");
            }
            break;
        case 3:
            if (mpi_rank == 0) {
                printf("The system is a MHD + test particle system.\n");
            }
            break;
        default:
            if (mpi_rank == 0) {
                printf("Error system option: [%d]\n", system_type);
            }
            return -1;
    }
    free(buff);
    return system_type;
}

/******************************************************************************
 * Read special information for the used system. 
 *  1. Loop current and wire current information for the wire-loop
 *     current system.
 *  2. Velocity field and force-free field information for force-free field
 *     system.
 *  3. Velocity field and magnetic field from MHD simulation for MHD +
 *     test-particle simulation. 
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  system_type: the system type.
 *  config_file_name: the configuration file name.
 ******************************************************************************/
int get_system_info(int mpi_rank, int system_type, char *config_file_name)
{
    int multi_tframe;
    switch (system_type) {
        case 0:
            if (mpi_rank == 0) {
                printf("Test system. All fields are calculated analytically.");
            }
            break;
        case 1:
            /* Wire-loop current system. */
            read_wlcs(mpi_rank, config_file_name);
            break;
        /* case 2: */
        /*     double v0; */
        /*     int tot_grid; */
        /*     dim_vfield(&v0); */
        /*     get_param_ff(); */
        /*     tot_grid = simul_grid.nx*simul_grid.ny*simul_grid.nz; */
        /*     vfd_b = (struct vfields*)malloc(sizeof(struct vfields)*tot_grid); */
        /*     read_vfields(0, v0, vfd_b); */
        /*     if (multi_tframe == 1) { */
        /*         vfd_a = (struct vfields*)malloc(sizeof(struct vfields)*tot_grid); */
        /*         read_vfields(1, v0, vfd_a); */
        /*     } */
        default:
            if (mpi_rank == 0) {
                printf("Error system option: [%d]\n", system_type);
            }
            return -1;
    }
    return 0;
}
