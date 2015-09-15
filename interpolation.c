#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "constants.h"
#include "domain.h"

struct grids simul_grid;
struct domain simul_domain;
int multi_tframe;

/******************************************************************************
 * set shared variables for interpolation procedures.
 ******************************************************************************/
void set_variabls_interpolation(grids *sgrid, domain *sdomain, int mt)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    multi_tframe = mt;
}

/******************************************************************************
 * Get the grid indices for the actual positions x, y, z, t. It will return two
 * sets of indices which contains the point (x, y, z, t) and the displacements
 * from those grid indices.
 *
 * Input:
 *  x, y, z, t: spatial positions and time.
 *
 * Output:
 *  ix1, iy1, iz1, it1: left grid indices.
 *  ix2, iy2, iz2, it2: right grid indices.
 *  dx, dy, dz, dt: displacements from the left grid indices.
 ******************************************************************************/
void grid_indices(double x, double y, double z, double t, int *ix1, int *iy1, 
        int *iz1, int *it1, int *ix2, int *iy2, int *iz2, int *it2,
        double *dx, double *dy, double *dz, double *dt)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    //nt = simul_grid.nt;
    /* left and right grid indices */
    *ix1 = floor(x/simul_grid.dx);
    *iy1 = floor(y/simul_grid.dy);
    *iz1 = floor(z/simul_grid.dz);
    *ix2 = *ix1 + 1;
    *iy2 = *iy1 + 1;
    *iz2 = *iz1 + 1;
    *dx = x/simul_grid.dx - (*ix1);
    *dy = y/simul_grid.dy - (*iy1);
    *dz = z/simul_grid.dz - (*iz1);
   
    /* boundary conditions */
    if (x > simul_domain.xmax) {
        *ix1 = nx - 1;
        *ix2 = 0;
        *dx = x/simul_grid.dx - (nx-1);
    }
    if (y > simul_domain.ymax) {
        *iy1 = ny - 1;
        *iy2 = 0;
        *dy = y/simul_grid.dy - (ny-1);
    }
    if (z > simul_domain.zmax) {
        *iz1 = nz - 1;
        *iz2 = 0;
        *dz = z/simul_grid.dz - (nz-1);
    }
    if (x < simul_domain.xmin) {
        *ix1 = nx - 1;
        *ix2 = 0;
        *dx = x/simul_grid.dx + 1.0;
    }
    if (y < simul_domain.ymin) {
        *iy1 = ny - 1;
        *iy2 = 0;
        *dy = y/simul_grid.dy + 1.0;
    }
    if (z < simul_domain.zmin) {
        *iz1 = nz - 1;
        *iz2 = 0;
        *dz = z/simul_grid.dz + 1.0;
    }

    if (multi_tframe == 0) {
        *it1 = 0;
        *it2 = 0;
        *dt = 0.0;
    }
    else if (multi_tframe == 1) {
        *it1 = floor(t/simul_grid.dt);
        *it2 = *it1 + 1;
        *dt = t/simul_grid.dt - (*it1);
    } else {
        printf("ERROR: wrong choice of time slices.\n");
        exit(1);
    }
}

/******************************************************************************
 * Transfer the indices of multidimensional array to the index on a one
 * dimensional array.
 *
 * Input:
 *  dims: the sizes along each dimension.
 *  indices: the grid indices to transformed to 1D.
 *  rank: the rank of the grid.
 *
 * Output:
 *  index: the 1D index.
 ******************************************************************************/
void index_transfer(int *dims, int *indices, int rank, long int *index)
{
    int i, shift;
    *index = indices[rank-1];
    shift = dims[rank-1];
    for (i = rank - 2; i >= 0; i--) {
        *index += shift * indices[i];
        shift *= dims[i];
    }
}

/******************************************************************************
 * Get the weights and corners indicies for trilinear interpolation.
 *
 * Input:
 *  x, y, z, t: spatial positions and time.
 *
 * Output:
 *  corners: the 1D indices for the corners.
 *  weights: the weights for trilinear interpolation.
 *  dt: the shift from the previous time point.
 ******************************************************************************/
void get_trilinar_parameters(double x, double y, double z, double t,
        long int *corners, double *weights, double *dt)
{
    int ix1, iy1, iz1, it1;
    int ix2, iy2, iz2, it2;
    double dx, dy, dz;
    int dims[3], indices[3];

    grid_indices(x, y, z, t, &ix1, &iy1, &iz1, &it1,
            &ix2, &iy2, &iz2, &it2, &dx, &dy, &dz, dt);

    dims[0] = simul_grid.nz; 
    dims[1] = simul_grid.ny; 
    dims[2] = simul_grid.nx;

    indices[0] = iz1; indices[1] = iy1; indices[2] = ix1;
    index_transfer(dims, indices, 3, &corners[0]);
    indices[0] = iz2; indices[1] = iy1; indices[2] = ix1;
    index_transfer(dims, indices, 3, &corners[1]);
    indices[0] = iz1; indices[1] = iy2; indices[2] = ix1;
    index_transfer(dims, indices, 3, &corners[2]);
    indices[0] = iz1; indices[1] = iy1; indices[2] = ix2;
    index_transfer(dims, indices, 3, &corners[3]);
    indices[0] = iz1; indices[1] = iy2; indices[2] = ix2;
    index_transfer(dims, indices, 3, &corners[4]);
    indices[0] = iz2; indices[1] = iy1; indices[2] = ix2;
    index_transfer(dims, indices, 3, &corners[5]);
    indices[0] = iz2; indices[1] = iy2; indices[2] = ix1;
    index_transfer(dims, indices, 3, &corners[6]);
    indices[0] = iz2; indices[1] = iy2; indices[2] = ix2;
    index_transfer(dims, indices, 3, &corners[7]);

    weights[0] = (1.0-dx) * (1.0-dy) * (1.0-dz);
    weights[1] = dx * (1.0-dy) * (1.0-dz);
    weights[2] = (1.0-dx) * dy * (1.0-dz);
    weights[3] = (1.0-dx) * (1.0-dy) * dz;
    weights[4] = (1.0-dx) * dy * dz;
    weights[5] = dx * (1.0-dy) * dz;
    weights[6] = dx * dy * (1.0-dz);
    weights[7] = dx * dy * dz;
}
