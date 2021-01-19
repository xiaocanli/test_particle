/******************************************************************************
* This file is part of CHAOTICB.
* Copyright (C) <2012-2014> <Xiaocan Li> <xl0009@uah.edu>
* 
* CHAOTICB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published ey
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
#include "electric_field.h"
#include "interpolation.h"

struct grids simul_grid;
struct domain simul_domain;
// 'a' is ahead, 'b' is behind, 'c' is current.
struct efields_double *efa_double, *efb_double, *efc_double;
struct efields_float *efa_float, *efb_float, *efc_float;
// sizeof(float) for float, sizeof(double) for double
int data_type;
double E0;
int multi_tframe;

/******************************************************************************
 * Set the global variable for current file.
 *
 * Input:
 *  sgrid: the grid information for the fields.
 *  sdomain: the domain information for the simulation.
 *  E0_field: the normalization for the electric field.
 *  mt: flag for whether multiple time slices of the fields are used.
 *
 * Variables to set:
 *  simul_grid: the grid information for the fields.
 *  simul_domain: the domain information for the simulation.
 *  E0: the normalization for the velocity field.
 ******************************************************************************/
void set_variables_efield(grids *sgrid, domain *sdomain, double E0_field,
        int mt, int dtype)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    E0 = E0_field;
    multi_tframe = mt;
    data_type = dtype;
}

/******************************************************************************
 * Initialize the electric field. The fields data are in double.
 ******************************************************************************/
void initialize_efield_double(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        efa_double = (efields_double *)malloc(sizeof(efields_double)*nx*ny*nz);
        efb_double = (efields_double *)malloc(sizeof(efields_double)*nx*ny*nz);
    } else {
        efc_double = (efields_double *)malloc(sizeof(efields_double)*nx*ny*nz);
    }
}

/******************************************************************************
 * Initialize the electric field. The fields data are in float.
 ******************************************************************************/
void initialize_efield_float(void)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    if (multi_tframe == 1) {
        efa_float = (efields_float *)malloc(sizeof(efields_float)*nx*ny*nz);
        efb_float = (efields_float *)malloc(sizeof(efields_float)*nx*ny*nz);
    } else {
        efc_float = (efields_float *)malloc(sizeof(efields_float)*nx*ny*nz);
    }
}

/******************************************************************************
 * Initialize electric field.
 ******************************************************************************/
void initialize_efield(void)
{
    if (data_type == sizeof(float)) {
        initialize_efield_float();
    }
    else if (data_type == sizeof(double)) {
        initialize_efield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Free the electric field with double data.
 ******************************************************************************/
void free_efield_double(void)
{
    if (multi_tframe == 1) {
        free(efa_double);
        free(efb_double);
    } else {
        free(efc_double);
    }
}

/******************************************************************************
 * Free the electric field with float data.
 ******************************************************************************/
void free_efield_float(void)
{
    if (multi_tframe == 1) {
        free(efa_float);
        free(efb_float);
    } else {
        free(efc_float);
    }
}

/******************************************************************************
 * Free the electric field.
 ******************************************************************************/
void free_efield(void)
{
    if (data_type == sizeof(float)) {
        free_efield_float();
    }
    else if (data_type == sizeof(double)) {
        free_efield_double();
    } else {
        printf("ERROR: wrong data type.\n");
        exit(1);
    }
}

/******************************************************************************
 * Read the electric fields from binary file. The data is in double.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_efields_double_binary(char *filepath, int ct)
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
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_double[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_double[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", ct, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_double[i].ez = data[i] * E0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_double[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_double[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_double[i].ez = data[i] * E0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_double[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_double[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", tindex, ".gda");
        read_data_serial_double_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_double[i].ez = data[i] * E0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the electric fields from binary file. The data is in float.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 ******************************************************************************/
void read_efields_float_binary(char *filepath, int ct)
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
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_float[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_float[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", ct, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efc_float[i].ez = data[i] * E0;
        }
    } else {
        // previous time point
        if (ct > 0) {
            tindex = ct - 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_float[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_float[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efb_float[i].ez = data[i] * E0;
        }

        // next time point.
        if (ct < simul_grid.nt - 1) {
            tindex = ct + 1;
        } else {
            tindex = ct;
        }
        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ex_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_float[i].ex = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ey_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_float[i].ey = data[i] * E0;
        }

        snprintf(fname, sizeof(fname), "%s%s%d%s", filepath, "ez_", tindex, ".gda");
        read_data_serial_float_binary(fname, offset, size, data);
        for (i = 0; i < nx*ny*nz; i++) {
            efa_float[i].ez = data[i] * E0;
        }
    }

    free(data);
}

/******************************************************************************
 * Read the electric fields from binary file.
 *
 * Input:
 *  filepath: the file path includes the files.
 *  ct: the current time slice.
 *  data_type: float or double.
 ******************************************************************************/
void read_efields_binary(char *filepath, int ct, int data_type)
{
    if (data_type == sizeof(float)) {
        read_efields_float_binary(filepath, ct);
    } else if (data_type == sizeof(double)) {
        read_efields_double_binary(filepath, ct);
    } else {
        printf("ERROR: wrong data type\n");
        exit(1);
    }
}

/******************************************************************************
 * Get electric field at a point. The fields have double data type.
 *
 * Input:
 *  x, y, z: the spatial positions.
 *  t: current time.
 ******************************************************************************/
void get_double_efield_at_point(double x, double y, double z, double t,
        double *ex, double *ey, double *ez)
{
    long int corners[8];
    double weights[8];
    double dt, ex1, ex2, ey1, ey2, ez1, ez2;
    *ex = 0.0; *ey = 0.0; *ez = 0.0;
    ex1 = 0.0; ey1 = 0.0; ez1 = 0.0;
    ex2 = 0.0; ey2 = 0.0; ez2 = 0.0;
    get_trilinar_parameters(x, y, z, t, corners, weights, &dt);
    if (multi_tframe == 0) {
        for (int i = 0; i < 8; i++) {
            *ex += efc_double[corners[i]].ex * weights[i];
            *ey += efc_double[corners[i]].ey * weights[i];
            *ez += efc_double[corners[i]].ez * weights[i];
        }
    } else {
        for (int i = 0; i < 8; i++) {
            ex1 += efa_double[corners[i]].ex * weights[i];
            ey1 += efa_double[corners[i]].ey * weights[i];
            ez1 += efa_double[corners[i]].ez * weights[i];
            ex2 += efb_double[corners[i]].ex * weights[i];
            ey2 += efb_double[corners[i]].ey * weights[i];
            ez2 += efb_double[corners[i]].ez * weights[i];
        }
        *ex = dt * ex1 + (1.0 - dt) * ex2;
        *ey = dt * ey1 + (1.0 - dt) * ey2;
        *ez = dt * ez1 + (1.0 - dt) * ez2;
    }
}

/******************************************************************************
 * Get electric field at a point. The fields have float data type.
 *
 * Input:
 *  x, y, z: the spatial positions.
 *  t: current time.
 ******************************************************************************/
void get_float_efield_at_point(double x, double y, double z, double t,
        float *ex, float *ey, float *ez)
{
    long int corners[8];
    double weights[8];
    double dt, ex1, ex2, ey1, ey2, ez1, ez2;
    *ex = 0.0; *ey = 0.0; *ez = 0.0;
    ex1 = 0.0; ey1 = 0.0; ez1 = 0.0;
    ex2 = 0.0; ey2 = 0.0; ez2 = 0.0;
    get_trilinar_parameters(x, y, z, t, corners, weights, &dt);
    if (multi_tframe == 0) {
        for (int i = 0; i < 8; i++) {
            *ex += efc_float[corners[i]].ex * weights[i];
            *ey += efc_float[corners[i]].ey * weights[i];
            *ez += efc_float[corners[i]].ez * weights[i];
        }
    } else {
        for (int i = 0; i < 8; i++) {
            ex1 += efa_float[corners[i]].ex * weights[i];
            ey1 += efa_float[corners[i]].ey * weights[i];
            ez1 += efa_float[corners[i]].ez * weights[i];
            ex2 += efb_float[corners[i]].ex * weights[i];
            ey2 += efb_float[corners[i]].ey * weights[i];
            ez2 += efb_float[corners[i]].ez * weights[i];
        }
        *ex = dt * ex1 + (1.0 - dt) * ex2;
        *ey = dt * ey1 + (1.0 - dt) * ey2;
        *ez = dt * ez1 + (1.0 - dt) * ez2;
    }
}
