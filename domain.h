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

int read_domain(int mpi_rank, char *config_file_name, domain *simul_domain,
        int *system_type, int *tracking_method);
