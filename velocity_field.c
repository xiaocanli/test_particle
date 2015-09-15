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
#include "domain.h"
#include "data_io.h"
#include "velocity_field.h"

struct grids simul_grid;
struct domain simul_domain;
struct vfields *vfd; 
double v0;

/******************************************************************************
 * Set the global variable for current file.
 *
 * Input:
 *  sgrid: the grid information for the fields.
 *  sdomain: the domain information for the simulation.
 *  v0_field: the normalization for the velocity field.
 *
 * Variables to set:
 *  simul_grid: the grid information for the fields.
 *  simul_domain: the domain information for the simulation.
 *  v0: the normalization for the velocity field.
 ******************************************************************************/
void set_variables_velocity(grids *sgrid, domain *sdomain, double v0_field)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    v0 = v0_field;
}

/******************************************************************************
 * Initialize the velocity field.
 ******************************************************************************/
void initialize_vfield(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    vfd = (vfields *)malloc(nx * ny * nz * sizeof(vfields));
}

/******************************************************************************
 * Free the velocity field.
 ******************************************************************************/
void free_vfield(void)
{
    free(vfd);
}

/******************************************************************************
 * Read the velocity fields from HDF5 files.
 *
 * Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 ******************************************************************************/
void read_vfields_h5(int ct, char *fname, char *gname)
{
    char dset_name[16];
    int rank = 4;
    int nz, ny, nx;
    long int i;

    hsize_t count[rank], offset[rank];

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    double *data = (double *)malloc(sizeof(double)*nz*ny*nx);

    count[0] = 1; count[1] = nz;
    count[2] = ny; count[3] = nx;
    offset[0] = ct; offset[1] = 0;
    offset[2] = 0; offset[3] = 0;
    snprintf(dset_name, sizeof(dset_name), "%s", "ux4d");
    read_data_serial_h5(rank, count, offset, fname, gname, dset_name, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vx = data[i] * v0;
    }
    snprintf(dset_name, sizeof(dset_name), "%s", "uy4d");
    read_data_serial_h5(rank, count, offset, fname, gname, dset_name, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vy = data[i] * v0;
    }
    snprintf(dset_name, sizeof(dset_name), "%s", "uz4d");
    read_data_serial_h5(rank, count, offset, fname, gname, dset_name, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vz = data[i] * v0;
    }
    free(data);
}

/******************************************************************************
 * Read the velocity fields from binary file.
 *
 * Input:
 *  ct: the current time slice.
 ******************************************************************************/
void read_vfields_binary(int ct)
{
    char fname[LEN_MAX];
    int nz, ny, nx;
    long int i, offset, size;

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    float *data = (float *)malloc(sizeof(float)*nz*ny*nx);

    offset = 0;
    size = nx * ny * nz;

    snprintf(fname, sizeof(fname), "%s%d%s", "vx_", ct, ".gda");
    read_data_serial_binary(fname, offset, size, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vx = data[i] * v0;
    }

    snprintf(fname, sizeof(fname), "%s%d%s", "vy_", ct, ".gda");
    read_data_serial_binary(fname, offset, size, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vy = data[i] * v0;
    }

    snprintf(fname, sizeof(fname), "%s%d%s", "vz_", ct, ".gda");
    read_data_serial_binary(fname, offset, size, data);
    for (i = 0; i < nx*ny*nz; i++) {
        vfd[i].vz = data[i] * v0;
    }

    free(data);
}
