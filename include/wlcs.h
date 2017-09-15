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

#ifndef WLCS_H
#define WLCS_H

/* Wire-loop current systems */
typedef struct wlcs{
	double x_wr, y_wr, z_wr;
	double costheta_wr, cosphi_wr;
	double sintheta_wr, sinphi_wr;
    double cur_wr, t0;
	double x_lp, y_lp, z_lp, r_lp;
	double cosalpha_lp, cosbeta_lp;
	double sinalpha_lp, sinbeta_lp;
    double cur_lp;
    double omega;
} wlcs;

#endif

#ifndef I0
#define I0 4.0E8
#endif

void read_wlcs(int mpi_rank, char *config_file_name);
void getemf_wlcs(double x, double y, double z, double t, 
        struct emfields *emf_tot);
void free_config();
void adjust_dt_normI(double *dt);
