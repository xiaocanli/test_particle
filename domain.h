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

void read_domain(int my_id, domain *simul_domain);
