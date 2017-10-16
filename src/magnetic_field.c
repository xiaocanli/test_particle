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
#include "magnetic_field.h"
#include "interpolation.h"

struct grids simul_grid;
struct domain simul_domain;
// 'a' is ahead, 'b' is behind, 'c' is current.
struct bfields_double *maga_double, *magb_double, *magc_double;
struct bfields_float *maga_float, *magb_float, *magc_float;
// sizeof(float) for float, sizeof(double) for double
int data_type;
double B0;
int multi_tframe;

/******************************************************************************
 * Set the global variable for current file.
 *
 * Input:
 *  sgrid: the grid information for the fields.
 *  sdomain: the domain information for the simulation.
 *  B0_field: the normalization for the magnetic field.
 *  mt: flag for whether multiple time slices of the fields are used.
 *
 * Variables to set:
 *  simul_grid: the grid information for the fields.
 *  simul_domain: the domain information for the simulation.
 *  B0: the normalization for the magnetic field.
 ******************************************************************************/
void set_variables_bfield(grids *sgrid, domain *sdomain, double B0_field,
        int mt, int dtype)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    B0 = B0_field;
    multi_tframe = mt;
    data_type = dtype;
}

/******************************************************************************
 * Initialize the magnetic field. The fields data are in double.
 ******************************************************************************/
void initialize_bfield_double(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        maga_double = (bfields_double *)malloc(nx*ny*nz * sizeof(bfields_double));
        magb_double = (bfields_double *)malloc(nx*ny*nz * sizeof(bfields_double));
    } else {
        magc_double = (bfields_double *)malloc(nx*ny*nz * sizeof(bfields_double));
    }
}

/******************************************************************************
 * Initialize the magnetic field. The fields data are in float.
 ******************************************************************************/
void initialize_bfield_float(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        maga_float = (bfields_float *)malloc(nx*ny*nz * sizeof(bfields_float));
        magb_float = (bfields_float *)malloc(nx*ny*nz * sizeof(bfields_float));
    } else {
        magc_float = (bfields_float *)malloc(nx*ny*nz * sizeof(bfields_float));
    }
}

/******************************************************************************
 * Initialize magnetic field.
 ******************************************************************************/
void initialize_bfield(void)
{
    if (data_type == sizeof(float)) {
        initialize_bfield_float();
    }
    else if (data_type == sizeof(double)) {
        initialize_bfield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Free the magnetic field with double data.
 ******************************************************************************/
void free_bfield_double(void)
{
    if (multi_tframe == 1) {
        free(maga_double);
        free(magb_double);
    } else {
        free(magc_double);
    }
}

/******************************************************************************
 * Free the magnetic field with float data.
 ******************************************************************************/
void free_bfield_float(void)
{
    if (multi_tframe == 1) {
        free(maga_float);
        free(magb_float);
    } else {
        free(magc_float);
    }
}

/******************************************************************************
 * Free the magnetic field.
 ******************************************************************************/
void free_bfield(void)
{
    if (data_type == sizeof(float)) {
        free_bfield_float();
    }
    else if (data_type == sizeof(double)) {
        free_bfield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Read the magnetic fields from HDF5 files. The data is in double.
 *
 o Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 ******************************************************************************/
void read_bfields_double_h5(int ct, char *fname, char *gname)
{
    char dset_bx[16], dset_by[16], dset_bz[16];
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

    snprintf(dset_bx, sizeof(dset_bx), "%s", "bx4d");
    snprintf(dset_by, sizeof(dset_by), "%s", "by4d");
    snprintf(dset_bz, sizeof(dset_bz), "%s", "bz4d");

    if (multi_tframe == 0) {
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].bx = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].by = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].bz = data[i] * B0;
        }
    } else {
        // fields at previous time point.
        if (ct > 0) {
            offset[0] = ct - 1;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].bx = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].by = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].bz = data[i] * B0;
        }

        // fields at next time point.
        if (ct < simul_grid.nt - 1) {
            offset[0] = ct + 1;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].bx = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].by = data[i] * B0;
        }
        read_data_serial_double_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].bz = data[i] * B0;
        }
    }
    free(data);
}

/******************************************************************************
 * Read the magnetic fields from HDF5 files. The data is in float.
 *
 o Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 ******************************************************************************/
void read_bfields_float_h5(int ct, char *fname, char *gname)
{
    char dset_bx[16], dset_by[16], dset_bz[16];
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

    snprintf(dset_bx, sizeof(dset_bx), "%s", "bx4d");
    snprintf(dset_by, sizeof(dset_by), "%s", "by4d");
    snprintf(dset_bz, sizeof(dset_bz), "%s", "bz4d");

    if (multi_tframe == 0) {
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].bx = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].by = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].bz = data[i] * B0;
        }
    } else {
        // fields at previous time point.
        if (ct > 0) {
            offset[0] = ct - 1;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].bx = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].by = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].bz = data[i] * B0;
        }

        // fields at next time point.
        if (ct < simul_grid.nt - 1) {
            offset[0] = ct + 1;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bx, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].bx = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_by, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].by = data[i] * B0;
        }
        read_data_serial_float_h5(rank, count, offset, fname, gname, dset_bz, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].bz = data[i] * B0;
        }
    }
    free(data);
}

/******************************************************************************
 * Read the magnetic fields from HDF5 files.
 *
 * Input:
 *  ct: time slice of the data.
 *  fname: the HDF5 file name.
 *  gname: the group name.
 *  data_type: the data type (float or double).
 ******************************************************************************/
void read_bfields_h5(int ct, char *fname, char *gname, int data_type)
{
    if (data_type == sizeof(float)) {
        read_bfields_float_h5(ct, fname, gname);
    } else if (data_type == sizeof(double)) {
        read_bfields_double_h5(ct, fname, gname);
    } else {
        printf("ERROR: wrong data type\n");
        exit(1);
    }
}

/******************************************************************************
 * Read the magnetic fields from binary file. The data is in double.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_bfields_double_binary(char *filepath, int ct)
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
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_double[i].bz = data[i] * B0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_double[i].bz = data[i] * B0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_double[i].bz = data[i] * B0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the magnetic fields from binary file. The data is in float.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_bfields_float_binary(char *filepath, int ct)
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
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magc_float[i].bz = data[i] * B0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            magb_float[i].bz = data[i] * B0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bx_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].bx = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "by_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].by = data[i] * B0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "bz_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            maga_float[i].bz = data[i] * B0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the magnetic fields from binary file.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 *  data_type: float or double.
 ******************************************************************************/
void read_bfields_binary(char *filepath, int ct, int data_type)
{
    if (data_type == sizeof(float)) {
        read_bfields_float_binary(filepath, ct);
    } else if (data_type == sizeof(double)) {
        read_bfields_double_binary(filepath, ct);
    } else {
        printf("ERROR: wrong data type\n");
        exit(1);
    }
}

/******************************************************************************
 * Get magnetic field at a point. The fields have double data type.
 *
 * Input:
 *  x, y, z: the spatial positions.
 *  t: current time.
 ******************************************************************************/
void get_double_bfield_at_point(double x, double y, double z, double t,
        double *bx, double *by, double *bz)
{
    long int corners[8];
    double weights[8];
    double dt, bx1, bx2, by1, by2, bz1, bz2;
    *bx = 0.0; *by = 0.0; *bz = 0.0;
    bx1 = 0.0; by1 = 0.0; bz1 = 0.0;
    bx2 = 0.0; by2 = 0.0; bz2 = 0.0;
    get_trilinar_parameters(x, y, z, t, corners, weights, &dt);
    if (multi_tframe == 0) {
        for (int i = 0; i < 8; i++) {
            *bx += magc_double[corners[i]].bx * weights[i];
            *by += magc_double[corners[i]].by * weights[i];
            *bz += magc_double[corners[i]].bz * weights[i];
        }
    } else {
        for (int i = 0; i < 8; i++) {
            bx1 += maga_double[corners[i]].bx * weights[i];
            by1 += maga_double[corners[i]].by * weights[i];
            bz1 += maga_double[corners[i]].bz * weights[i];
            bx2 += magb_double[corners[i]].bx * weights[i];
            by2 += magb_double[corners[i]].by * weights[i];
            bz2 += magb_double[corners[i]].bz * weights[i];
        }
        *bx = dt * bx1 + (1.0 - dt) * bx2;
        *by = dt * by1 + (1.0 - dt) * by2;
        *bz = dt * bz1 + (1.0 - dt) * bz2;
    }
}

/******************************************************************************
 * Get magnetic field at a point. The fields have float data type.
 *
 * Input:
 *  x, y, z: the spatial positions.
 *  t: current time.
 ******************************************************************************/
void get_float_bfield_at_point(double x, double y, double z, double t,
        float *bx, float *by, float *bz)
{
    long int corners[8];
    double weights[8];
    double dt, bx1, bx2, by1, by2, bz1, bz2;
    *bx = 0.0; *by = 0.0; *bz = 0.0;
    bx1 = 0.0; by1 = 0.0; bz1 = 0.0;
    bx2 = 0.0; by2 = 0.0; bz2 = 0.0;
    get_trilinar_parameters(x, y, z, t, corners, weights, &dt);
    if (multi_tframe == 0) {
        for (int i = 0; i < 8; i++) {
            *bx += magc_float[corners[i]].bx * weights[i];
            *by += magc_float[corners[i]].by * weights[i];
            *bz += magc_float[corners[i]].bz * weights[i];
        }
    } else {
        for (int i = 0; i < 8; i++) {
            bx1 += maga_float[corners[i]].bx * weights[i];
            by1 += maga_float[corners[i]].by * weights[i];
            bz1 += maga_float[corners[i]].bz * weights[i];
            bx2 += magb_float[corners[i]].bx * weights[i];
            by2 += magb_float[corners[i]].by * weights[i];
            bz2 += magb_float[corners[i]].bz * weights[i];
        }
        *bx = dt * bx1 + (1.0 - dt) * bx2;
        *by = dt * by1 + (1.0 - dt) * by2;
        *bz = dt * bz1 + (1.0 - dt) * bz2;
    }
}
