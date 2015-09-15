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
// 'a' is ahead, 'b' is behind, 'c' is current.
struct vfields_double *va_double, *vb_double, *vc_double;
struct vfields_float *va_float, *vb_float, *vc_float;
// sizeof(float) for float, sizeof(double) for double
int data_type;
double v0;
int multi_tframe;

/******************************************************************************
 * Set the global variable for current file.
 *
 * Input:
 *  sgrid: the grid information for the fields.
 *  sdomain: the domain information for the simulation.
 *  v0_field: the normalization for the velocity field.
 *  mt: flag for whether multiple time slices of the fields are used.
 *
 * Variables to set:
 *  simul_grid: the grid information for the fields.
 *  simul_domain: the domain information for the simulation.
 *  v0: the normalization for the velocity field.
 ******************************************************************************/
void set_variables_velocity(grids *sgrid, domain *sdomain, double v0_field,
        int mt, int dtype)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    v0 = v0_field;
    multi_tframe = mt;
    data_type = dtype;
}

/******************************************************************************
 * Initialize the velocity field. The fields data are in double.
 ******************************************************************************/
void initialize_vfield_double(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        va_double = (vfields_double *)malloc(nx*ny*nz * sizeof(vfields_double));
        vb_double = (vfields_double *)malloc(nx*ny*nz * sizeof(vfields_double));
    } else {
        vc_double = (vfields_double *)malloc(nx*ny*nz * sizeof(vfields_double));
    }
}

/******************************************************************************
 * Initialize the velocity field. The fields data are in float.
 ******************************************************************************/
void initialize_vfield_float(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        va_float = (vfields_float *)malloc(nx*ny*nz * sizeof(vfields_float));
        vb_float = (vfields_float *)malloc(nx*ny*nz * sizeof(vfields_float));
    } else {
        vc_float = (vfields_float *)malloc(nx*ny*nz * sizeof(vfields_float));
    }
}

/******************************************************************************
 * Initialize velocity field.
 ******************************************************************************/
void initialize_vfield(void)
{
    if (data_type == sizeof(float)) {
        initialize_vfield_float();
    }
    else if (data_type == sizeof(double)) {
        initialize_vfield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Free the velocity field with double data.
 ******************************************************************************/
void free_vfield_double(void)
{
    if (multi_tframe == 1) {
        free(va_double);
        free(vb_double);
    } else {
        free(vc_double);
    }
}

/******************************************************************************
 * Free the velocity field with float data.
 ******************************************************************************/
void free_vfield_float(void)
{
    if (multi_tframe == 1) {
        free(va_float);
        free(vb_float);
    } else {
        free(vc_float);
    }
}

/******************************************************************************
 * Free the velocity field.
 ******************************************************************************/
void free_vfield(void)
{
    if (data_type == sizeof(float)) {
        free_vfield_float();
    }
    else if (data_type == sizeof(double)) {
        free_vfield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Read the velocity fields from HDF5 files. The data is in double.
 *
 * Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 ******************************************************************************/
void read_vfields_double_h5(int ct, char *fname, char *gname)
{
    char dset_ux[16], dset_uy[16], dset_uz[16];
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

    snprintf(dset_ux, sizeof(dset_ux), "%s", "ux4d");
    snprintf(dset_uy, sizeof(dset_uy), "%s", "uy4d");
    snprintf(dset_uz, sizeof(dset_uz), "%s", "uz4d");

    if (multi_tframe == 0) {
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vx = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vy = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vz = data[i] * v0;
        }
    } else {
        // fields at previous time point.
        if (ct > 0) {
            offset[0] = ct - 1;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vx = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vy = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vz = data[i] * v0;
        }

        // fields at next time point.
        if (ct < simul_grid.nt - 1) {
            offset[0] = ct + 1;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vx = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vy = data[i] * v0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vz = data[i] * v0;
        }
    }
    free(data);
}

/******************************************************************************
 * Read the velocity fields from HDF5 files. The data is in float.
 *
 * Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 ******************************************************************************/
void read_vfields_float_h5(int ct, char *fname, char *gname)
{
    char dset_ux[16], dset_uy[16], dset_uz[16];
    int rank = 4;
    int nz, ny, nx;
    long int i;

    hsize_t count[rank], offset[rank];

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    float *data = (float *)malloc(sizeof(float)*nz*ny*nx);

    count[0] = 1; count[1] = nz;
    count[2] = ny; count[3] = nx;
    offset[0] = ct; offset[1] = 0;
    offset[2] = 0; offset[3] = 0;

    snprintf(dset_ux, sizeof(dset_ux), "%s", "ux4d");
    snprintf(dset_uy, sizeof(dset_uy), "%s", "uy4d");
    snprintf(dset_uz, sizeof(dset_uz), "%s", "uz4d");

    if (multi_tframe == 0) {
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vx = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vy = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vz = data[i] * v0;
        }
    } else {
        // fields at previous time point.
        if (ct > 0) {
            offset[0] = ct - 1;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vx = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vy = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vz = data[i] * v0;
        }

        // fields at next time point.
        if (ct < simul_grid.nt - 1) {
            offset[0] = ct + 1;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_ux, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vx = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uy, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vy = data[i] * v0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_uz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vz = data[i] * v0;
        }
    }
    free(data);
}

/******************************************************************************
 * Read the velocity fields from HDF5 files.
 *
 * Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 *  data_type: the data type (float or double).
 ******************************************************************************/
void read_vfields_h5(int ct, char *fname, char *gname, int data_type)
{
    if (data_type == sizeof(float)) {
        read_vfields_float_h5(ct, fname, gname);
    } else if (data_type == sizeof(double)) {
        read_vfields_double_h5(ct, fname, gname);
    } else {
        printf("ERROR: wrong data type\n");
        exit(1);
    }
}

/******************************************************************************
 * Read the velocity fields from binary file. The data is in double.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_vfields_double_binary(char *filepath, int ct)
{
    char fname[LEN_MAX];
    int nz, ny, nx;
    long int i, offset, size;
    int tindex;

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    double *data = (double *)malloc(sizeof(double)*nz*ny*nx);

    offset = 0;
    size = nx * ny * nz;

    if (multi_tframe == 0) {
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_double[i].vz = data[i] * v0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_double[i].vz = data[i] * v0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_double[i].vz = data[i] * v0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the velocity fields from binary file. The data is in float.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_vfields_float_binary(char *filepath, int ct)
{
    char fname[LEN_MAX];
    int nz, ny, nx;
    long int i, offset, size;
    int tindex;

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    float *data = (float *)malloc(sizeof(float)*nz*ny*nx);

    offset = 0;
    size = nx * ny * nz;

    if (multi_tframe == 0) {
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vc_float[i].vz = data[i] * v0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            vb_float[i].vz = data[i] * v0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vx_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vx = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vy_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vy = data[i] * v0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "vz_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            va_float[i].vz = data[i] * v0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the velocity fields from binary file.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 *  data_type: float or double.
 ******************************************************************************/
void read_vfields_binary(char *filepath, int ct, int data_type)
{
    if (data_type == sizeof(float)) {
        read_vfields_float_binary(filepath, ct);
    } else if (data_type == sizeof(double)) {
        read_vfields_double_binary(filepath, ct);
    } else {
        printf("ERROR: wrong data type\n");
        exit(1);
    }
}
