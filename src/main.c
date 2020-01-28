#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include "constants.h"
#include "diagnostics.h"
#include "wlcs.h"
#include "electric_field.h"
#include "magnetic_field.h"
#include "poincare_map.h"
#include "test_particle.h"
#include "velocity_field.h"

int get_run_type(int mpi_rank, char *config_file_name, int *run_type);

/******************************************************************************
 * Main program for test particle simulation.
 ******************************************************************************/
int main(int argc, char **argv)
{
    int ierr, mpi_size, mpi_rank;
    char config_file_name[LEN_MAX];
    int ipvd, err, system_type, tracking_method, run_type;

    /* ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ipvd); */
    /* ierr = MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &ipvd); */
    ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &ipvd);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    snprintf(config_file_name, sizeof(config_file_name), "%s", "init.dat");
    get_run_type(mpi_rank, config_file_name, &run_type);
    struct domain simul_domain;
    err = read_domain(mpi_rank, config_file_name, run_type, &simul_domain,
            &system_type, &tracking_method);
    if (err < 0) {
        ierr = MPI_Finalize();
        return 0;
    }

    if (run_type == 1) {
        test_particle(mpi_rank, mpi_size, config_file_name, system_type,
                tracking_method, simul_domain);
    }
    else if (run_type == 2) {
        poincare_map(mpi_rank, mpi_size, config_file_name, &simul_domain);
    }

    switch (system_type) {
        case 0:
            break;
        case 1:
            free_config();
            break;
        case 2:
            free_vfield();
            break;
        case 3:
            free_vfield();
            free_bfield();
            break;
        case 4:
            if (run_type == 1)
                free_efield();
            free_bfield();
            break;
        default:
            if (mpi_rank == 0) {
                printf("Error system option: [%d]\n", system_type);
            }
            ierr = MPI_Finalize();
            return -1;
    }
    ierr = MPI_Finalize();
    return 0;
}

/******************************************************************************
 * Get the run type from the configuration file
 ******************************************************************************/
int get_run_type(int mpi_rank, char *config_file_name, int *run_type)
{
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    int err, msg;
    FILE *fp;

    fp = fopen(config_file_name, "r");

    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Simulation run type") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "Run type (1: test particle simulation 2: Poincare map):%d\n", run_type);
    if (msg != 1) {
        printf("Failed to read the simulation run type.\n");
        exit(1);
    }
    switch (*run_type) {
        case 1:
            if (mpi_rank == 0) {
                printf("This is a test-particle simulation.\n");
            }
            break;
        case 2:
            if (mpi_rank == 0) {
                printf("This simulation gets a Poincare map.\n");
            }
            break;
        default:
            if (mpi_rank == 0) {
                printf("ERROR: wrong flag for run type.");
            }
            return -1;
    }
    fclose(fp);

    free(buff);
    return 0;
}
