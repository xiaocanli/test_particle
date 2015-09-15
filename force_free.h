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
#include "emfields.h"

#ifndef ABC_FIELD_H
#define ABC_FIELD_H
#define ffA 1.0
#define ffB sqrt(2.0/3.0)
#define ffC sqrt(1.0/3.0)
#endif

/* Parameters for force-free field */
#ifndef PARAM_FF_H
#define PARAM_FF_H
typedef struct param_ff {
    double B0, lambda0;
} param_ff;
#endif

/* Velocity fields */
#ifndef VFIELDS_H
#define VFIELDS_H
typedef struct vfields {
	double vx, vy, vz;
} vfields;
#endif

void getb_ff(double x, double y, double z, double t, struct bfields *bmf);

void gete_ff(double x, double y, double z, double t, struct bfields *bmf,
        struct efields *elf);

void getemf_ff(double x, double y, double z, double t, struct emfields *emf);

void get_param_ff(int mpi_rank, char *config_file_name);

void set_variables_ff(grids *sgrid, domain *sdomain);
