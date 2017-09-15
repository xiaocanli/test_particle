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
#include "domain.h"

void read_ptlinfo(int *nptl_tot, double *vthe);
void gaussian_rand(double *ran1, double *ran2);
void init_ptl(int nptl, double vthe, domain simul_domain, struct particles *ptl);
void getemf_test(double x, double y, double z, double t, struct emfields *emf);
void get_emf(double x, double y, double z, double t, struct emfields *emf);
void assign_ptl(int *nptl, double *vthe);
void release_memory(struct particles *ptl);
void get_emf(double x, double y, double z, double t, struct emfields *emf);
void tracking_time_method(double *tott);
