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
#include "velocity_field.h"

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
    getb_ff(x, y, z, t, &bmf);
    gete_ff(x, y, z, t, &bmf, &elf);
    emf->Bx = bmf.Bx;
    emf->By = bmf.By;
    emf->Bz = bmf.Bz;
    emf->Ex = elf.Ex;
    emf->Ey = elf.Ey;
    emf->Ez = elf.Ez;
}

/******************************************************************************
 * Calculate the inductive electric fields from turbulent velocity fields.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 *  bmf: the magnetic field.
 *
 * Output:
 *  elf: 3 components of electric field.
 ******************************************************************************/
void gete_ff(double x, double y, double z, double t, struct bfields *bmf,
        struct efields *elf)
{
    double vx, vy, vz;
    get_double_velocity_at_point(x, y, z, t, &vx, &vy, &vz);
    elf->Ex = bmf->By*vz - bmf->Bz*vy;
    elf->Ey = bmf->Bz*vx - bmf->Bx*vz;
    elf->Ez = bmf->Bx*vy - bmf->By*vx;
}

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
