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
#include "Global.h"
#include "force_free.h"

int imultiple;
struct grids simul_grid;
struct param_ff param;

void ReadDataSerialHDF5(int rank, hsize_t *count, 
        hsize_t *offset, char *fname, char *gname, char *dset_name, 
        double *data);

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
void grid_indices(double x, double y, double z, double t, int *ix1, int *iy1, 
        int *iz1, int *it1, int *ix2, int *iy2, int *iz2, int *it2,
        double *dx, double *dy, double *dz, double *dt)
{
    int nx, ny, nz;
    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    //nt = simul_grid.nt;
    /* left and right grid indices */
    *ix1 = floor(x/simul_grid.dx);
    *iy1 = floor(y/simul_grid.dy);
    *iz1 = floor(z/simul_grid.dz);
    *ix2 = *ix1 + 1;
    *iy2 = *iy1 + 1;
    *iz2 = *iz1 + 1;
    *dx = x/simul_grid.dx - (*ix1);
    *dy = y/simul_grid.dy - (*iy1);
    *dz = z/simul_grid.dz - (*iz1);
   
    /* boundary conditions */
    if (x > simul_domain.xmax) {
        *ix1 = nx - 1;
        *ix2 = 0;
        *dx = x/simul_grid.dx - (nx-1);
    }
    if (y > simul_domain.ymax) {
        *iy1 = ny - 1;
        *iy2 = 0;
        *dy = y/simul_grid.dy - (ny-1);
    }
    if (z > simul_domain.zmax) {
        *iz1 = nz - 1;
        *iz2 = 0;
        *dz = z/simul_grid.dz - (nz-1);
    }
    if (x < simul_domain.xmin) {
        *ix1 = nx - 1;
        *ix2 = 0;
        *dx = x/simul_grid.dx + 1.0;
    }
    if (y < simul_domain.ymin) {
        *iy1 = ny - 1;
        *iy2 = 0;
        *dy = y/simul_grid.dy + 1.0;
    }
    if (z < simul_domain.zmin) {
        *iz1 = nz - 1;
        *iz2 = 0;
        *dz = z/simul_grid.dz + 1.0;
    }

    if (imultiple == 0) {
        *it1 = 0;
        *it2 = 0;
        *dt = 0.0;
    }
    else if (imultiple == 1) {
        *it1 = floor(t/simul_grid.dt);
        *it2 = *it1 + 1;
        *dt = t/simul_grid.dt - (*it1);
    }
}

/******************************************************************************
 * Calculate the inductive electric fields from turbulent velocity fields.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 * Output:
 *  elf: 3 components of electric field.
 ******************************************************************************/
void gete_ff(double x, double y, double z, double t, 
        struct bfields *bmf, struct efields *elf)
{
    int ix1, iy1, iz1, it1;
    int ix2, iy2, iz2, it2;
    double dx, dy, dz, dt;
    int i000, i100, i010, i001, i101, i011, i110, i111;
    double v1, v2, v3, v4, v5, v6, v7, v8;
    int dims[3], indices[3];
    double vx, vy, vz;
    double vx1, vy1, vz1;
    grid_indices(x, y, z, t, &ix1, &iy1, &iz1, &it1,
            &ix2, &iy2, &iz2, &it2, &dx, &dy, &dz, &dt);
    /* Trilinear  interpolation */
    dims[0] = simul_grid.nz; 
    dims[1] = simul_grid.ny; 
    dims[2] = simul_grid.nx;
    indices[0] = iz1; indices[1] = iy1; indices[2] = ix1;
    index_transfer(dims, indices, 3, &i000);
    indices[0] = iz2; indices[1] = iy1; indices[2] = ix1;
    index_transfer(dims, indices, 3, &i001);
    indices[0] = iz1; indices[1] = iy2; indices[2] = ix1;
    index_transfer(dims, indices, 3, &i010);
    indices[0] = iz1; indices[1] = iy1; indices[2] = ix2;
    index_transfer(dims, indices, 3, &i100);
    indices[0] = iz1; indices[1] = iy2; indices[2] = ix2;
    index_transfer(dims, indices, 3, &i011);
    indices[0] = iz2; indices[1] = iy1; indices[2] = ix2;
    index_transfer(dims, indices, 3, &i101);
    indices[0] = iz2; indices[1] = iy2; indices[2] = ix1;
    index_transfer(dims, indices, 3, &i110);
    indices[0] = iz2; indices[1] = iy2; indices[2] = ix2;
    index_transfer(dims, indices, 3, &i111);

    v1 = (1.0-dx)*(1.0-dy)*(1.0-dz);
    v2 = dx*(1.0-dy)*(1.0-dz);
    v3 = (1.0-dx)*dy*(1.0-dz);
    v4 = (1.0-dx)*(1.0-dy)*dz;
    v5 = (1.0-dx)*dy*dz;
    v6 = dx*(1.0-dy)*dz;
    v7 = dx*dy*(1.0-dz);
    v8 = dx*dy*dz;

    vx = vfd_b[i000].vx*v1 + vfd_b[i001].vx*v2 + vfd_b[i010].vx*v3 +
         vfd_b[i100].vx*v4 + vfd_b[i011].vx*v5 + vfd_b[i101].vx*v6 +
         vfd_b[i110].vx*v7 + vfd_b[i111].vx*v8;
    vy = vfd_b[i000].vy*v1 + vfd_b[i001].vy*v2 + vfd_b[i010].vy*v3 +
         vfd_b[i100].vy*v4 + vfd_b[i011].vy*v5 + vfd_b[i101].vy*v6 +
         vfd_b[i110].vy*v7 + vfd_b[i111].vy*v8;
    vz = vfd_b[i000].vz*v1 + vfd_b[i001].vz*v2 + vfd_b[i010].vz*v3 +
         vfd_b[i100].vz*v4 + vfd_b[i011].vz*v5 + vfd_b[i101].vz*v6 +
         vfd_b[i110].vz*v7 + vfd_b[i111].vz*v8;

    if (imultiple == 1) {
        /* Interpolate along t, but can do spatial interpolation first,
         * since it is linear interpolation.
         */
        vx1 = vfd_a[i000].vx*v1 + vfd_a[i001].vx*v2 + vfd_a[i010].vx*v3 +
             vfd_a[i100].vx*v4 + vfd_a[i011].vx*v5 + vfd_a[i101].vx*v6 +
             vfd_a[i110].vx*v7 + vfd_a[i111].vx*v8;
        vy1 = vfd_a[i000].vy*v1 + vfd_a[i001].vy*v2 + vfd_a[i010].vy*v3 +
             vfd_a[i100].vy*v4 + vfd_a[i011].vy*v5 + vfd_a[i101].vy*v6 +
             vfd_a[i110].vy*v7 + vfd_a[i111].vy*v8;
        vz1 = vfd_a[i000].vz*v1 + vfd_a[i001].vz*v2 + vfd_a[i010].vz*v3 +
             vfd_a[i100].vz*v4 + vfd_a[i011].vz*v5 + vfd_a[i101].vz*v6 +
             vfd_a[i110].vz*v7 + vfd_a[i111].vz*v8;

        vx = (1.0-dt)*vx + dt*vx1;
        vy = (1.0-dt)*vy + dt*vy1;
        vz = (1.0-dt)*vz + dt*vz1;
    }

    elf->Ex = bmf->By*vz - bmf->Bz*vy;
    elf->Ey = bmf->Bz*vx - bmf->Bx*vz;
    elf->Ez = bmf->Bx*vy - bmf->By*vx;
}

/******************************************************************************
 * Transfer the indices of multidimensional array to the index on a one
 * dimensional array.
 ******************************************************************************/
void index_transfer(int *dims, int *indices, int ndim, int *index)
{
    int i, shift;
    *index = indices[ndim-1];
    shift = dims[ndim-1];
    for (i=ndim-2; i>=0; i--) {
        *index += shift*indices[i];
        shift *= dims[i];
    }
}

/******************************************************************************
 * Read the normalized B and the inverse of the length scale of the magnetic
 * fields from the initialization file.
 *
 * Output:
 *  param: the parameters mentioned above.
 ******************************************************************************/
void get_param_ff(void)
{
    FILE *fp;
    char buff[LEN_MAX];
    int msg;
    fp = fopen("init.dat", "r");
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
    if (my_id == 0) {
        printf("The normalized magnetic field is %lf Gauss\n", param.B0);
        printf("1/Length scale (lambda) is %lf \n", param.lambda0);
    }
}

/******************************************************************************
 * Calculate the magnetic fields from analytical expressions for force-free
 * magnetic field.
 *
 * Input:
 *  x, y, z, t: current spatial positions and time.
 * Output:
 *  bmf: 3 components of magnetic field.
 ******************************************************************************/
void getb_ff(double x, double y, double z, double t, struct bfields *bmf)
{
    double B0, lambda0;
    B0 = param.B0;
    lambda0 = param.lambda0;
    /* One solution of \nabla\times\vec{B} = \lambda\vect{B} 
     * lamada0 here is constant.
     */
    bmf->Bx = ffA*sin(lambda0*z) + ffC*cos(lambda0*y);
    bmf->By = ffB*sin(lambda0*x) + ffA*cos(lambda0*z);
    bmf->Bz = ffC*sin(lambda0*y) + ffB*cos(lambda0*x);

    bmf->Bx *= B0;
    bmf->By *= B0;
    bmf->Bz *= B0;
}

/******************************************************************************
 * Read the dimensions of velocity fields generated from inverse FFT.
 *
 * Output:
 *  simul_grid: structure contains grid dimensions and grid sizes.
 ******************************************************************************/
void dim_vfield(double *v0)
{
    FILE *fp;
    char buff[LEN_MAX];
    int msg;
    fp = fopen("init.dat", "r");
    while (fgets(buff,LEN_MAX,fp) != NULL){
        if (strstr(buff, "Velocity fields")) {
            break;
        }
    };
    msg = fscanf(fp, "x dimension(nx): %d\n", &simul_grid.nx);
    if (msg != 1) {
        printf("Failed to read nx.\n");
        exit(1);
    }
    msg = fscanf(fp, "y dimension(ny): %d\n", &simul_grid.ny);
    if (msg != 1) {
        printf("Failed to read ny.\n");
        exit(1);
    }
    msg = fscanf(fp, "z dimension(nz): %d\n", &simul_grid.nz);
    if (msg != 1) {
        printf("Failed to read nz.\n");
        exit(1);
    }
    msg = fscanf(fp, "time slices(nt): %d\n", &simul_grid.nt);
    if (msg != 1) {
        printf("Failed to read nt.\n");
        exit(1);
    }
    msg = fscanf(fp, "Normalized v0:%lf\n", v0);
    if (msg != 1) {
        printf("Failed to read normalized v0.\n");
        exit(1);
    }
    msg = fscanf(fp, "If multiple t slice: %d\n", &imultiple);
    if (msg != 1) {
        printf("Failed to read the flag of whether to have multiple"
                "time slices.\n");
        exit(1);
    }
    fclose(fp);
    /* Grid sizes */
    simul_grid.dx = (simul_domain.xmax-simul_domain.xmin)/(simul_grid.nx-1.0);
    simul_grid.dy = (simul_domain.ymax-simul_domain.ymin)/(simul_grid.ny-1.0);
    simul_grid.dz = (simul_domain.zmax-simul_domain.zmin)/(simul_grid.nz-1.0);
    simul_grid.dt = (simul_domain.tmax-simul_domain.tmin)/(simul_grid.nt-1.0);
    if (my_id == 0) {
        printf("Dimensions of velocity fields (nx, ny, nz, nt): "
                "(%d, %d, %d, %d)\n", simul_grid.nx, 
                simul_grid.ny, simul_grid.nz, simul_grid.nt );
        printf("Normalized velocity: %lf\n", *v0);
    }
}

/******************************************************************************
 * Read the velocity fields generated from inverse FFT from HDF5 files.
 *
 * Input:
 *  nx, ny, nz: dimensions of the data, which are global data.
 *  it: time slice of the data.
 * Output:
 *  vfd: structure array of the velocity fields.
 ******************************************************************************/
void read_vfields(int it, double v0, struct vfields *vfd)
{
    char *fname, *gname, dset_name[256];
    int rank = 4;
    int nz, ny, nx;
    long int i;

    hsize_t count[rank], offset[rank];

    nx = simul_grid.nx;
    ny = simul_grid.ny;
    nz = simul_grid.nz;
    double *data = (double *)malloc(sizeof(double)*nz*ny*nx);

    fname = "data/u4d.h5";
    gname = "/u4d";
    count[0] = 1; count[1] = nz;
    count[2] = ny; count[3] = nx;
    offset[0] = it; offset[1] = 0;
    offset[2] = 0; offset[3] = 0;
    snprintf(dset_name, sizeof(dset_name), "%s", "ux4d");
    ReadDataSerialHDF5(rank, count, offset, fname, gname,
            dset_name, data);
    for (i=0; i<nx*ny*nz; i++) {
        vfd[i].vx = data[i] * v0;
    }
    snprintf(dset_name, sizeof(dset_name), "%s", "uy4d");
    ReadDataSerialHDF5(rank, count, offset, fname, gname,
            dset_name, data);
    for (i=0; i<nx*ny*nz; i++) {
        vfd[i].vy = data[i] * v0;
    }
    snprintf(dset_name, sizeof(dset_name), "%s", "uz4d");
    ReadDataSerialHDF5(rank, count, offset, fname, gname,
            dset_name, data);
    for (i=0; i<nx*ny*nz; i++) {
        vfd[i].vz = data[i] * v0;
    }
    free(data);
}

/******************************************************************************
 * Read data from HDF5 file.
 *
 * Input:
 *  count: the local dimensions of the data in current MPI process.
 *  offset: the local offset in each dimension in current MPI process.
 *  fname: the name of the HDF5 file.
 *  gname: the name of the group.
 *  dset_name: dataset name.
 * Output:
 *  data: the read data from the file.
 ******************************************************************************/
void ReadDataSerialHDF5(int rank, hsize_t *count, hsize_t *offset,
        char *fname, char *gname, char *dset_name, double *data)
{
    hid_t file_id, group_id, dset_id;
    hid_t filespace, memspace;
    //herr_t status;

    /* Open the existing file. */
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Open a group */
    group_id = H5Gopen1(file_id, gname);

    /* Open a dataset and get its dataspace */
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    filespace = H5Dget_space(dset_id);

    memspace = H5Screate_simple(rank, count, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, 
            offset, NULL, count, NULL);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace,
            filespace, H5P_DEFAULT, data);

    /* Close/release resources */
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Gclose(group_id);
    H5Fclose(file_id);
}
