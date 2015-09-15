#include "domain.h"

void set_variabls_interpolation(grids *sgrid, domain *sdomain, int mt);

void grid_indices(double x, double y, double z, double t, int *ix1, int *iy1, 
        int *iz1, int *it1, int *ix2, int *iy2, int *iz2, int *it2,
        double *dx, double *dy, double *dz, double *dt);

void get_trilinar_parameters(double x, double y, double z, double t,
        long int *corners, double *weights, double *dt);
