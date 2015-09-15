#ifndef DOMAIN_H
#define DOMAIN_H

/* Simulation domain */
typedef struct domain{
    double xmin, xmax;
    double ymin, ymax;
    double zmin, zmax;
    double tmin, tmax;
} domain;

#endif

/* Domain information for grids */
#ifndef GRIDS_H
#define GRIDS_H

typedef struct grids{
    int nx, ny, nz, nt;
    double dx, dy, dz, dt;
} grids;

#endif

int read_domain(int mpi_rank, char *config_file_name, domain *simul_domain,
        int *system_type, int *tracking_method);

void get_fields_dims(int mpi_rank, char *config_file_name, domain simul_domain,
        grids *simul_grid, double *v0_field, double *B0_field, int *multi_tframe);
