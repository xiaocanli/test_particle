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
#define ffA 1.0
#define ffB sqrt(2.0/3.0)
#define ffC sqrt(1.0/3.0)

void getb_ff(double x, double y, double z, double t, struct bfields *bmf);
void gete_ff(double x, double y, double z, double t, 
        struct bfields *bmf, struct efields *elf);
void grid_indices(double x, double y, double z, double t, int *ix1, int *iy1, 
        int *iz1, int *it1, int *ix2, int *iy2, int *iz2, int *it2,
        double *dx, double *dy, double *dz, double *dt);
void index_transfer(int *dims, int *indices, int ndim, int *index);
void getemf_ff(double x, double y, double z, double t, struct emfields *emf);
void dim_vfield(double *v0);
void read_vfields(int it, double v0, struct vfields *vfd);
void get_param_ff(void);
