/******************************************************************************
* This file is part of CHAOTICB.
* Copyright (C) <2012-2014> <Xiaocan Li> <xl0009@uah.edu>
* 
* CHAOTICB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* CHAOTICB is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with CHAOTICB.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <hdf5.h>
#include "constants.h"
#include "force_free.h"
#include "emfields.h"
#include "domain.h"

struct grids simul_grid;
struct domain simul_domain;
struct param_ff param;

void ReadDataSerialHDF5(int rank, hsize_t *count, hsize_t *offset, char *fname,
        char *gname, char *dset_name, double *data);

/******************************************************************************
 * Get the electromagnetic fields for a force-free field case. The magnetic
 * field is calculated from analytical expressions. The electric field is
 * calculated from velocity field by calculating the inductive E field.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 * Output:
 *  emf: electromagnetic fields.
 ******************************************************************************/
void getemf_ff(double x, double y, double z, double t, struct emfields *emf)
{
    struct bfields bmf;
    struct efields elf;
    /* getb_ff(x, y, z, t, &bmf); */
    /* gete_ff(x, y, z, t, &bmf, &elf); */
    emf->Bx = bmf.Bx;
    emf->By = bmf.By;
    emf->Bz = bmf.Bz;
    emf->Ex = elf.Ex;
    emf->Ey = elf.Ey;
    emf->Ez = elf.Ez;
}

/******************************************************************************
 * Get the grid indices for the actual positions x, y, z, t. It will return
 * two sets of indices which contains the point (x, y, z, t) and the
 * displacements from those grid indices.
 *
 * Input:
 *  x, y, z, t: spatial positions and time.
 * Output:
 *  ix1, iy1, iz1, it1: left grid indices.
 *  ix2, iy2, iz2, it2: right grid indices.
 *  dx, dy, dz, dt: displacements from the left grid indices.
 ******************************************************************************/
/* void grid_indices(double x, double y, double z, double t, int *ix1, int *iy1, */ 
/*         int *iz1, int *it1, int *ix2, int *iy2, int *iz2, int *it2, */
/*         double *dx, double *dy, double *dz, double *dt) */
/* { */
/*     int nx, ny, nz; */
/*     nx = simul_grid.nx; */
/*     ny = simul_grid.ny; */
/*     nz = simul_grid.nz; */
/*     //nt = simul_grid.nt; */
/*     /1* left and right grid indices *1/ */
/*     *ix1 = floor(x/simul_grid.dx); */
/*     *iy1 = floor(y/simul_grid.dy); */
/*     *iz1 = floor(z/simul_grid.dz); */
/*     *ix2 = *ix1 + 1; */
/*     *iy2 = *iy1 + 1; */
/*     *iz2 = *iz1 + 1; */
/*     *dx = x/simul_grid.dx - (*ix1); */
/*     *dy = y/simul_grid.dy - (*iy1); */
/*     *dz = z/simul_grid.dz - (*iz1); */
   
/*     /1* boundary conditions *1/ */
/*     if (x > simul_domain.xmax) { */
/*         *ix1 = nx - 1; */
/*         *ix2 = 0; */
/*         *dx = x/simul_grid.dx - (nx-1); */
/*     } */
/*     if (y > simul_domain.ymax) { */
/*         *iy1 = ny - 1; */
/*         *iy2 = 0; */
/*         *dy = y/simul_grid.dy - (ny-1); */
/*     } */
/*     if (z > simul_domain.zmax) { */
/*         *iz1 = nz - 1; */
/*         *iz2 = 0; */
/*         *dz = z/simul_grid.dz - (nz-1); */
/*     } */
/*     if (x < simul_domain.xmin) { */
/*         *ix1 = nx - 1; */
/*         *ix2 = 0; */
/*         *dx = x/simul_grid.dx + 1.0; */
/*     } */
/*     if (y < simul_domain.ymin) { */
/*         *iy1 = ny - 1; */
/*         *iy2 = 0; */
/*         *dy = y/simul_grid.dy + 1.0; */
/*     } */
/*     if (z < simul_domain.zmin) { */
/*         *iz1 = nz - 1; */
/*         *iz2 = 0; */
/*         *dz = z/simul_grid.dz + 1.0; */
/*     } */

/*     if (imultiple == 0) { */
/*         *it1 = 0; */
/*         *it2 = 0; */
/*         *dt = 0.0; */
/*     } */
/*     else if (imultiple == 1) { */
/*         *it1 = floor(t/simul_grid.dt); */
/*         *it2 = *it1 + 1; */
/*         *dt = t/simul_grid.dt - (*it1); */
/*     } */
/* } */

/******************************************************************************
 * Calculate the inductive electric fields from turbulent velocity fields.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 * Output:
 *  elf: 3 components of electric field.
 ******************************************************************************/
/* void gete_ff(double x, double y, double z, double t, */ 
/*         struct bfields *bmf, struct efields *elf) */
/* { */
/*     int ix1, iy1, iz1, it1; */
/*     int ix2, iy2, iz2, it2; */
/*     double dx, dy, dz, dt; */
/*     int i000, i100, i010, i001, i101, i011, i110, i111; */
/*     double v1, v2, v3, v4, v5, v6, v7, v8; */
/*     int dims[3], indices[3]; */
/*     double vx, vy, vz; */
/*     double vx1, vy1, vz1; */
/*     grid_indices(x, y, z, t, &ix1, &iy1, &iz1, &it1, */
/*             &ix2, &iy2, &iz2, &it2, &dx, &dy, &dz, &dt); */
/*     /1* Trilinear  interpolation *1/ */
/*     dims[0] = simul_grid.nz; */ 
/*     dims[1] = simul_grid.ny; */ 
/*     dims[2] = simul_grid.nx; */
/*     indices[0] = iz1; indices[1] = iy1; indices[2] = ix1; */
/*     index_transfer(dims, indices, 3, &i000); */
/*     indices[0] = iz2; indices[1] = iy1; indices[2] = ix1; */
/*     index_transfer(dims, indices, 3, &i001); */
/*     indices[0] = iz1; indices[1] = iy2; indices[2] = ix1; */
/*     index_transfer(dims, indices, 3, &i010); */
/*     indices[0] = iz1; indices[1] = iy1; indices[2] = ix2; */
/*     index_transfer(dims, indices, 3, &i100); */
/*     indices[0] = iz1; indices[1] = iy2; indices[2] = ix2; */
/*     index_transfer(dims, indices, 3, &i011); */
/*     indices[0] = iz2; indices[1] = iy1; indices[2] = ix2; */
/*     index_transfer(dims, indices, 3, &i101); */
/*     indices[0] = iz2; indices[1] = iy2; indices[2] = ix1; */
/*     index_transfer(dims, indices, 3, &i110); */
/*     indices[0] = iz2; indices[1] = iy2; indices[2] = ix2; */
/*     index_transfer(dims, indices, 3, &i111); */

/*     v1 = (1.0-dx)*(1.0-dy)*(1.0-dz); */
/*     v2 = dx*(1.0-dy)*(1.0-dz); */
/*     v3 = (1.0-dx)*dy*(1.0-dz); */
/*     v4 = (1.0-dx)*(1.0-dy)*dz; */
/*     v5 = (1.0-dx)*dy*dz; */
/*     v6 = dx*(1.0-dy)*dz; */
/*     v7 = dx*dy*(1.0-dz); */
/*     v8 = dx*dy*dz; */

/*     vx = vfd_b[i000].vx*v1 + vfd_b[i001].vx*v2 + vfd_b[i010].vx*v3 + */
/*          vfd_b[i100].vx*v4 + vfd_b[i011].vx*v5 + vfd_b[i101].vx*v6 + */
/*          vfd_b[i110].vx*v7 + vfd_b[i111].vx*v8; */
/*     vy = vfd_b[i000].vy*v1 + vfd_b[i001].vy*v2 + vfd_b[i010].vy*v3 + */
/*          vfd_b[i100].vy*v4 + vfd_b[i011].vy*v5 + vfd_b[i101].vy*v6 + */
/*          vfd_b[i110].vy*v7 + vfd_b[i111].vy*v8; */
/*     vz = vfd_b[i000].vz*v1 + vfd_b[i001].vz*v2 + vfd_b[i010].vz*v3 + */
/*          vfd_b[i100].vz*v4 + vfd_b[i011].vz*v5 + vfd_b[i101].vz*v6 + */
/*          vfd_b[i110].vz*v7 + vfd_b[i111].vz*v8; */

/*     if (imultiple == 1) { */
/*         /1* Interpolate along t, but can do spatial interpolation first, */
/*          * since it is linear interpolation. */
/*          *1/ */
/*         vx1 = vfd_a[i000].vx*v1 + vfd_a[i001].vx*v2 + vfd_a[i010].vx*v3 + */
/*              vfd_a[i100].vx*v4 + vfd_a[i011].vx*v5 + vfd_a[i101].vx*v6 + */
/*              vfd_a[i110].vx*v7 + vfd_a[i111].vx*v8; */
/*         vy1 = vfd_a[i000].vy*v1 + vfd_a[i001].vy*v2 + vfd_a[i010].vy*v3 + */
/*              vfd_a[i100].vy*v4 + vfd_a[i011].vy*v5 + vfd_a[i101].vy*v6 + */
/*              vfd_a[i110].vy*v7 + vfd_a[i111].vy*v8; */
/*         vz1 = vfd_a[i000].vz*v1 + vfd_a[i001].vz*v2 + vfd_a[i010].vz*v3 + */
/*              vfd_a[i100].vz*v4 + vfd_a[i011].vz*v5 + vfd_a[i101].vz*v6 + */
/*              vfd_a[i110].vz*v7 + vfd_a[i111].vz*v8; */

/*         vx = (1.0-dt)*vx + dt*vx1; */
/*         vy = (1.0-dt)*vy + dt*vy1; */
/*         vz = (1.0-dt)*vz + dt*vz1; */
/*     } */

/*     elf->Ex = bmf->By*vz - bmf->Bz*vy; */
/*     elf->Ey = bmf->Bz*vx - bmf->Bx*vz; */
/*     elf->Ez = bmf->Bx*vy - bmf->By*vx; */
/* } */

/******************************************************************************
 * Transfer the indices of multidimensional array to the index on a one
 * dimensional array.
 ******************************************************************************/
/* void index_transfer(int *dims, int *indices, int ndim, int *index) */
/* { */
/*     int i, shift; */
/*     *index = indices[ndim-1]; */
/*     shift = dims[ndim-1]; */
/*     for (i=ndim-2; i>=0; i--) { */
/*         *index += shift*indices[i]; */
/*         shift *= dims[i]; */
/*     } */
/* } */

/******************************************************************************
 * Read the normalized B and the inverse of the length scale of the magnetic
 * fields from the initialization file.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  config_file_name: the configuration filename.
 *
 * Output:
 *  param: the parameters for the force-free field.
 ******************************************************************************/
void get_param_ff(int mpi_rank, char *config_file_name)
{
    FILE *fp;
    char buff[LEN_MAX];
    int msg;
    fp = fopen(config_file_name, "r");
    while (fgets(buff,LEN_MAX,fp) != NULL){
        if (strstr(buff, "Force-free field")) {
            break;
        }
    };
    msg = fscanf(fp, "Normalized B(Gauss):%lf\n", &param.B0);
    if (msg != 1) {
        printf("Failed to read normalized B.\n");
        exit(1);
    }
    msg = fscanf(fp, "1/Length scale (lambda):%lf\n", &param.lambda0);
    if (msg != 1) {
        printf("Failed to read length scale lambda.\n");
        exit(1);
    }
    fclose(fp);
    if (mpi_rank == 0) {
        printf("The normalized magnetic field is %lf Gauss\n", param.B0);
        printf("1/Length scale (lambda) is %lf \n", param.lambda0);
    }
}

/******************************************************************************
 * Calculate the magnetic fields from analytical expressions for force-free
 * magnetic field. One solution of \nabla\times\vec{B} = \lambda\vect{B} 
 * Here, lamada0 is constant.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 *
 * Output:
 *  bmf: 3 components of magnetic field.
 ******************************************************************************/
void getb_ff(double x, double y, double z, double t, struct bfields *bmf)
{
    double B0, lambda0;
    B0 = param.B0;
    lambda0 = param.lambda0;
    bmf->Bx = ffA * sin(lambda0*z) + ffC * cos(lambda0*y);
    bmf->By = ffB * sin(lambda0*x) + ffA * cos(lambda0*z);
    bmf->Bz = ffC * sin(lambda0*y) + ffB * cos(lambda0*x);

    bmf->Bx *= B0;
    bmf->By *= B0;
    bmf->Bz *= B0;
}
