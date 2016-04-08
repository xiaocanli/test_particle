#include "domain.h"
#include "hdf5.h"

#ifndef PARTICLES_H
#define PARTICLES_H
typedef struct particles{
    double x, y, z;
    double vx, vy, vz, t;
} particles;
#endif

void assign_particles(int mpi_rank, int mpi_size, int *nptl_tot,
        int *nptl_traj_tot, int system_type, int *nptl, int *nptl_accumulate);

void read_particle_info(int mpi_rank, char *config_file_name, int *nptl_tot,
        double *ptemp, double *pmass, double *pcharge, int *charge_sign,
        double *vthe, double *charge_mass, int *bc_flag,
        int *is_traj_diagnostic, int *nptl_traj_tot, int *is_single_vel);

void particle_broadcast(int mpi_rank, int mpi_size, int nptl_tot, int *nptl,
        int *nptl_accumulate);

void initialize_partilces(int mpi_rank, int mpi_size, char *config_file_name,
        domain simul_domain, int nptl, double vthe, int system_type,
        int *nptl_accumulate,  particles *ptl, int is_single_vel);

void create_particles_ctype(hid_t *memtype, hid_t *filetype);

void read_particles(int mpi_rank, int nptl, int *nptl_accumulate,
        struct particles *ptl, char *fname);
