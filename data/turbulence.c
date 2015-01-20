#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "hdf5.h"
#include <mpi.h>
#include <fftw3-mpi.h>

#define MAX(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define M_PI 3.14159265358979323846
#define nquant 6
#define dimxyz 256
#define nslice 64

int ierr, num_procs, my_id;
int izs, ize, nzdim;

void frequency(int nx, int ny, int nz, int nt);
void index_transfer(int *dims, int *indices, int ndim, int *index);
void velocity(int nx, int ny, int nz, int nt, char *qname, double *data_max);
void normalization(int nx, int ny, int nz, int nt, double vmax, char *qname);
void InitHDF5Output(char *fname, char *gname, char *dnamelist[], int rank, 
        hsize_t *dimsf, hid_t *file_id, hid_t *group_id, hid_t *dset_id);
void WriteDataHDF5(hid_t dset_id, hid_t memspace, hid_t plist_id, 
        hsize_t *count, hsize_t *offset, double *data);
void GenerateRandomComplexVector(double kx, double ky, double kz, double l, 
        int ix, int iy, int iz, double *wx_re, double *wy_re, double *wz_re,
        double *wx_im, double *wy_im, double *wz_im);
void GenerateVelocityPhase(int nx, int ny, int nz, int nt, int it, double l, 
        double *vx_re, double *vy_re, double *vz_re, double *vx_im, 
        double *vy_im, double *vz_im);
void ReadDataParallelHDF5(int rank, hsize_t *dimsf, hsize_t *count, 
        hsize_t *offset, char *fname, char *gname, char *dset_name, 
        double *data);
void WriteDataParallelHDF5(int rank, hsize_t *dimsf, hsize_t *count, 
        hsize_t *offset, char *fname, char *gname, char *dset_name, 
        int newfile_flag,  int newgroup_flag, int newdataset_flag, double *data);
void GetLinearDistribution(long int sz, double dmin, double dmax,
        int nbands, char *qname, double *data);
void GetKineticEnergyDensity(int local_n0, int nt, int nz,
        int ny, int nx, double *data);

/******************************************************************************
 * Main function.
 ******************************************************************************/
int main(int argc, char **argv)
{
    int nx=dimxyz, ny=dimxyz; 
    int nz=dimxyz, nt=nslice;
    double vxmax, vymax, vzmax, vmax;

    ierr = MPI_Init(&argc, &argv);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &my_id);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    nzdim = nz / num_procs;
    izs = nzdim*my_id;
    ize = izs + nzdim;
    if (my_id < (nzdim % num_procs)) {
        ++nzdim;
        izs = izs + my_id;
        ize = izs + nzdim;
    }

    printf("id, izs, ize: %d %d %d\n", my_id, izs, ize);

    printf("Generating turbulent velocity field.\n");
    frequency(nx, ny, nz, nt);

    velocity(nx, ny, nz, nt, "ux", &vxmax);
    velocity(nx, ny, nz, nt, "uy", &vymax);
    velocity(nx, ny, nz, nt, "uz", &vzmax);
    vmax = MAX(vxmax, vymax);
    vmax = MAX(vmax, vzmax);

    /* Save the maximum velocity to file. */
    if (my_id==0) {
        FILE *fp;
        fp = fopen("vmax.dat", "w");
        fprintf(fp, "%lf %lf %lf %lf\n", vxmax, vymax, vzmax, vmax);
        fclose(fp);
    }

    normalization(nx, ny, nz, nt, vmax, "ux");
    normalization(nx, ny, nz, nt, vmax, "uy");
    normalization(nx, ny, nz, nt, vmax, "uz");

    ierr = MPI_Finalize();
    return 0;
}

/******************************************************************************
 * Generating turbulent velocity field in k-w space using the method from
 * 1. Stam, Jos, and Eugene Fiume. "Turbulent wind fields for gaseous phenomena." 
 * Proceedings of the 20th annual conference on Computer graphics and 
 * interactive techniques. ACM, 1993. APA  
 * 2. Stam, Jos. A general animation framework for gaseous phenomena. European 
 * Research Consortium for Informatics and Mathematics, 1997.
 *
 * Input:
 *  nx, ny, nz, nt: the dimensions of the data.
 ******************************************************************************/
void frequency(int nx, int ny, int nz, int nt)
{
    double l, kz;
    int it, iz;

    int rank=4;
    int iquant;
    char fname[] = "u4d_k.h5";
    char gname[] = "/u4d_k";
    char *dnamelist[nquant] = {"ux4d_re", "ux4d_im", "uy4d_re", 
        "uy4d_im", "uz4d_re", "uz4d_im"};
    hid_t file_id, group_id, dset_id[nquant];
    hid_t memspace;
    hid_t plist_id;
    //herr_t status;
    hsize_t dimsf[rank]; // dataset dimensions
    hsize_t count[rank], offset[rank];

    /* data space size for the dataset */
    dimsf[0] = nt;
    dimsf[1] = nz;
    dimsf[2] = ny;
    dimsf[3] = nx/2+1;

    /* Initialization for HDF5 output */
    InitHDF5Output(fname, gname, dnamelist, rank, 
            dimsf, &file_id, &group_id, dset_id);

    count[0] = 1;
    count[1] = nz;
    count[2] = ny;
    count[3] = nx/2+1;
    memspace = H5Screate_simple(rank, count, NULL);

    double *vx_re = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));
    double *vy_re = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));
    double *vz_re = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));
    double *vx_im = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));
    double *vy_im = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));
    double *vz_im = (double *) malloc(sizeof(double)*nz*ny*(nx/2+1));

    // Create property list for collective dataset write.
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    srand48((unsigned)time(NULL)+my_id*num_procs);
    //srand48(time(NULL));
    for (it=my_id; it<nt; it=it+num_procs) {
        printf("id, it: %d %d\n", my_id, it);
        l = ((double)it)/nt;
        if (l > nt/2) l -= 1.0;
        GenerateVelocityPhase(nx, ny, nz, nt, it, l, 
                vx_re, vy_re, vz_re, vx_im, vy_im, vz_im);
        offset[0] = it;
        offset[1] = 0;
        offset[2] = 0;
        offset[3] = 0;
        WriteDataHDF5(dset_id[0], memspace, plist_id, 
                count, offset, vx_re);
        WriteDataHDF5(dset_id[1], memspace, plist_id, 
                count, offset, vy_re);
        WriteDataHDF5(dset_id[2], memspace, plist_id, 
                count, offset, vz_re);
        WriteDataHDF5(dset_id[3], memspace, plist_id, 
                count, offset, vx_im);
        WriteDataHDF5(dset_id[4], memspace, plist_id, 
                count, offset, vy_im);
        WriteDataHDF5(dset_id[5], memspace, plist_id, 
                count, offset, vz_im);
    }

    free(vx_re); free(vy_re); free(vz_re);
    free(vx_im); free(vy_im); free(vz_im);

    /*
     * Close/release resources
     */
    for (iquant=0; iquant<nquant; iquant++) {
        H5Dclose(dset_id[iquant]);
    }
    //status = H5Sclose(filespace);
    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Procedure to write data to HDF5 file.
 * Input:
 *  dset_id, memspace, plist_id: handlers for dataset, memory space and
 *      property list.
 *  data: the actual data to write to the disk.
 ******************************************************************************/
void WriteDataHDF5(hid_t dset_id, hid_t memspace, hid_t plist_id, 
        hsize_t *count, hsize_t *offset, double *data)
{
    hid_t filespace;
    filespace = H5Dget_space(dset_id);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace,
            filespace, plist_id, data);
    H5Sclose(filespace);
}

/******************************************************************************
 * Generate random complex vector from Kolmogorov energy spectrum
 * and the temporal spread function.
 * 
 * Input:
 *  kx, ky, kz, l: wave numbers of temporal frequency.
 *  ix, iy, iz: current indices in the 3 dimensions.
 * Output:
 *  wx_re, wy_re, wz_re: real parts of the random complex vector.
 *  wx_im, wy_im, wz_im: imaginary parts of the random complex vector.
 ******************************************************************************/
void GenerateRandomComplexVector(double kx, double ky, double kz, double l, 
        int ix, int iy, int iz, double *wx_re, double *wy_re, double *wz_re,
        double *wx_im, double *wy_im, double *wz_im)
{
    double ks, k2, tx, ty, tz;
    double sfractal, ssmooth, skl, ramp;
    double twopi = 2.0*M_PI;
    k2 = kx*kx + ky*ky + kz*kz;
    ks = sqrt(k2);
    tx = drand48() * twopi;
    ty = drand48() * twopi;
    tz = drand48() * twopi;
    if ((ix==0) && (iy==0) && (iz==0)) {
        sfractal = 0.0; // Kolmogorov spectrum
        ssmooth = 0.0;  // Temporal spread function
    } else {
        sfractal = pow(ks, -5.0/3.0);
        ssmooth = exp(-l*l/k2)/ks;
        //ssmooth = exp(-l*l/k2);
        //ssmooth = 1.0;
    }
    skl = sfractal * ssmooth;
    ramp = skl;
    *wx_re = ramp * cos(tx);
    *wx_im = ramp * sin(tx);
    *wy_re = ramp * cos(ty);
    *wy_im = ramp * sin(ty);
    *wz_re = ramp * cos(tz);
    *wz_im = ramp * sin(tz);
}

/******************************************************************************
 * Generate velocity data in phase space for the energy spectrum using the
 * method from:
 *  1. Stam, Jos, and Eugene Fiume. "Turbulent wind fields for gaseous 
 *  phenomena." Proceedings of the 20th annual conference on Computer graphics 
 *  and interactive techniques. ACM, 1993. APA  
 *  2. Stam, Jos. A general animation framework for gaseous phenomena. European 
 *  Research Consortium for Informatics and Mathematics, 1997.
 * Here, a 2D data is generated for hdf5 output (x-y dimension).
 *
 * Input:
 *  nx, ny, nz: spatial dimensions of the data.
 *  iz, kz: the index in kz-dimension, normalized position iz/kz.
 *  l: normalized position in temporal frequency dimension.
 * Output:
 *  vx_re, vy_re, vz_re: real part of the velocity in phase space.
 *  vx_im, vy_im, vz_im: imaginary part of the velocity in phase space.
 ******************************************************************************/
void GenerateVelocityPhase(int nx, int ny, int nz, int nt, int it, double l, 
        double *vx_re, double *vy_re, double *vz_re, double *vx_im, 
        double *vy_im, double *vz_im)
{
    int ix, iy, iz;
    double kx, ky, kz, k2;
    double wx_re, wy_re, wz_re, wx_im, wy_im, wz_im;
    long int index_1d; /* 1D index for the 2D data for aligned memory */
    index_1d = 0;
    for (iz=0; iz<nz; iz++) {
        //printf("iz: %d\n", iz);
        kz = ((double)iz)/nz;
        if (iz > nz/2) kz -= 1.0;
        for (iy=0; iy<ny; iy++) {
            ky = ((double)iy)/ny;
            if (iy > ny/2) ky -= 1.0;
            for (ix=0; ix<=nx/2; ix++) {
                //printf("ix: %d\n", ix);
                kx = ((double)ix)/nx;
                k2 = kx*kx + ky*ky + kz*kz;
                GenerateRandomComplexVector(kx, ky, kz, l, 
                        ix, iy, iz, &wx_re, &wy_re, &wz_re,
                        &wx_im, &wy_im, &wz_im);

                // Project onto plane normal to {kx,ky,kz}
                if ((ix==0) && (iy==0) && (iz==0)) {
                    vx_re[index_1d] = 0.0;  
                    vy_re[index_1d] = 0.0; 
                    vz_re[index_1d] = 0.0; 
                    vx_im[index_1d] = 0.0;  
                    vy_im[index_1d] = 0.0; 
                    vz_im[index_1d] = 0.0; 
                } else {
                    vx_re[index_1d] = (1-kx*kx/k2)*wx_re - 
                        (kx*ky/k2)*wy_re - (kx*kz/k2)*wz_re; 
                    vy_re[index_1d] = (-ky*kx/k2)*wx_re + 
                        (1-ky*ky/k2)*wy_re - (ky*kz/k2)*wz_re;
                    vz_re[index_1d] = (-kz*kx/k2)*wx_re - 
                        (kz*ky/k2)*wy_re + (1-kz*kz/k2)*wz_re;
                    vx_im[index_1d] = (1-kx*kx/k2)*wx_im - 
                        (kx*ky/k2)*wy_im - (kx*kz/k2)*wz_im; 
                    vy_im[index_1d] = (-ky*kx/k2)*wx_im + 
                        (1-ky*ky/k2)*wy_im - (ky*kz/k2)*wz_im;
                    vz_im[index_1d] = (-kz*kx/k2)*wx_im - 
                        (kz*ky/k2)*wy_im + (1-kz*kz/k2)*wz_im;
                }

                // "mirror" symmetries at the axis of the symmetry
                if ((ix==0 || ix==nx/2) && (iy==0 || iy==ny/2) && 
                        (iz==0 || iz==nz/2) && (it==0 || it==nt/2)) {
                    vx_im[index_1d] = 0.0;
                    vy_im[index_1d] = 0.0;
                    vz_im[index_1d] = 0.0;
                }
                index_1d++;
            }
        }
    }
}

/******************************************************************************
 * Initialization for HDF5 file output.
 *
 * Input:
 *  fname, gname: file name and group name.
 *  dnamelist: dataset name list.
 *  rank: the dimensions of the data to save.
 *  dimsf: hdf5 datatype for integer arrays.
 * Output:
 *  file_id: the file handler for hdf5 file.
 *  group_id: the group handler for one group.
 *  dset_id: the dataset handlers for datasets. 
 ******************************************************************************/
void InitHDF5Output(char *fname, char *gname, char *dnamelist[], int rank, 
        hsize_t *dimsf, hid_t *file_id, hid_t *group_id, hid_t *dset_id)
{
    int iquant;
    hid_t filespace;
    hid_t plist_id;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // Disable ROMIO's data-sieving
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");
    // Enable ROMIO's collective buffering
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");

    // Setup file access property list with parallel I/O access.
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    MPI_Info_free(&info);
    // Create a new file collectively and release property list identifier.
    *file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    // Create a new group
    *group_id = H5Gcreate(*file_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // Create the dataspace for the dataset
    filespace = H5Screate_simple(rank, dimsf, NULL);
    // Create the dataset with default properties and close filespace.
    for (iquant=0; iquant<nquant; iquant++) {
        dset_id[iquant] = H5Dcreate(*group_id, dnamelist[iquant], 
                H5T_NATIVE_DOUBLE, filespace, H5P_DEFAULT, H5P_DEFAULT, 
                H5P_DEFAULT);
    }
    H5Sclose(filespace);
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
 * Read data from HDF5 file collectively.
 *
 * Input:
 *  dimsf: the dimension array of the data.
 *  count: the local dimensions of the data in current MPI process.
 *  offset: the local offset in each dimension in current MPI process.
 *  fname: the name of the HDF5 file.
 *  gname: the name of the group.
 *  dset_name: dataset name.
 * Output:
 *  data: the read data from the file.
 ******************************************************************************/
void ReadDataParallelHDF5(int rank, hsize_t *dimsf, hsize_t *count, 
        hsize_t *offset, char *fname, char *gname, char *dset_name, 
        double *data)
{
    hid_t file_id, group_id, dset_id;
    hid_t filespace, memspace;
    hid_t plist_id;
    //herr_t status;

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    // Disable ROMIO's data-sieving
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");
    // Enable ROMIO's collective buffering
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");

    // Setup file access property list with parallel I/O access.
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    MPI_Info_free(&info);
    // Open the existing file.
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, plist_id);
    H5Pclose(plist_id);
    // Open a group
    group_id = H5Gopen1(file_id, gname);

    // Open a dataset and get its dataspace
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    filespace = H5Dget_space(dset_id);

    // Create property list for collective dataset read/write
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    memspace = H5Screate_simple(rank, count, NULL);

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, data);

    /*
     * Close/release resources
     */
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Write data to HDF5 file collectively.
 *
 * Input:
 *  dimsf: the dimension array of the data.
 *  count: the local dimensions of the data in current MPI process.
 *  offset: the local offset in each dimension in current MPI process.
 *  fname: the name of the HDF5 file.
 *  gname: the name of the group.
 *  dset_name: dataset name.
 *  newfile_flag: the flag to check whether to create a new file.
 *  newgroup_flag: the flag to check whether to create a new group. 
 *  newdataset_flag: the flag to chekc whether to create a new dataset.
 *  data: the data to write to the hdf5 file.
 ******************************************************************************/
void WriteDataParallelHDF5(int rank, hsize_t *dimsf, hsize_t *count, 
        hsize_t *offset, char *fname, char *gname, char *dset_name, 
        int newfile_flag,  int newgroup_flag, int newdataset_flag, double *data)
{
    hid_t file_id, group_id, dset_id;
    hid_t filespace, memspace;
    hid_t plist_id;
    //herr_t status;

    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Info info = MPI_INFO_NULL;
    MPI_Info_create(&info);
    /* Disable ROMIO's data-sieving */
    MPI_Info_set(info, "romio_ds_read", "disable");
    MPI_Info_set(info, "romio_ds_write", "disable");
    /* Enable ROMIO's collective buffering */
    MPI_Info_set(info, "romio_cb_read", "enable");
    MPI_Info_set(info, "romio_cb_write", "enable");

    /* Setup file access property list with parallel I/O access. */
    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, comm, info);
    MPI_Info_free(&info);
    if (newfile_flag==1) {
        /* Create a new file collectively and release property list identifier. */
        file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    }
    else {
        /* Open an existing file. */
        file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    }
    H5Pclose(plist_id);
    if (newgroup_flag==1) {
        /* Create a new group */
        group_id = H5Gcreate(file_id, gname, 
                H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else {
        /* Open an existing group */
        group_id = H5Gopen1(file_id, gname);
    }

    if (newdataset_flag==1) {
        /* Create a dataspace for the dataset. */
        filespace = H5Screate_simple(rank, dimsf, NULL);
        /* Create a dataset with default properties. */
        dset_id = H5Dcreate(group_id, dset_name, H5T_NATIVE_DOUBLE,
                filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    }
    else {
        /* Open a dataset and get its dataspace */
        dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
        filespace = H5Dget_space(dset_id);
    }

    /* Create property list for collective dataset read/write */
    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    memspace = H5Screate_simple(rank, count, NULL);

    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, 
            filespace, plist_id, data);

    /*
     * Close/release resources
     */
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Pclose(plist_id);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * This procedure is for generating fields in real space.
 * Input:
 *  nx, ny, nt, nt: the dimensions of the data.
 *  qname: quantity name. It can be vx, vy, vz.
 * Output:
 *  data_max: the maximum of the absolute value of the data set.
 *  The generated data is written to the hard disk.
 ******************************************************************************/
void velocity(int nx, int ny, int nz, int nt, char *qname, double *data_max)
{
    int rank = 4;
    char buf1[256];
    double *data;
    char *fname, *gname, dset_name[256];
    ptrdiff_t len;
    const ptrdiff_t n[4] = {nt, nz, ny, nx};
    const ptrdiff_t nh[4] = {nt, nz, ny, nx/2+1};
    ptrdiff_t alloc_local, local_n0, local_0_start, i;
    double *rout;
    fftw_complex *cin;
    fftw_plan plan;

    fftw_mpi_init();
    /* get local data size and allocate */
    alloc_local = fftw_mpi_local_size(rank, nh, MPI_COMM_WORLD, 
            &local_n0, &local_0_start);
    cin = fftw_alloc_complex(alloc_local);
    rout = fftw_alloc_real(alloc_local*2);
    data = fftw_alloc_real(alloc_local);

    /* Create a plan for out-of-place c2r DFT */
    plan = fftw_mpi_plan_dft_c2r(rank, n, cin, rout, 
            MPI_COMM_WORLD, FFTW_ESTIMATE);

    /* read phase space velocity */
    hsize_t dimsf[rank], count[rank], offset[rank];
    fname = "u4d_k.h5";
    gname = "/u4d_k";
    dimsf[0] = nt; dimsf[1] = nz;
    dimsf[2] = ny; dimsf[3] = nx/2+1;
    count[0] = local_n0; count[1] = nz;
    count[2] = ny; count[3] = nx/2+1;
    offset[0] = local_0_start;
    offset[1] = 0; offset[2] = 0; offset[3] = 0;
    snprintf(dset_name, sizeof(dset_name), "%s%s", qname, "4d_re");
    ReadDataParallelHDF5(rank, dimsf, count, offset, 
            fname, gname, dset_name, data);
    len = local_n0*nh[1]*nh[2]*nh[3];
    for (i=0; i<len; i++) {
        cin[i][0] = data[i];
    }
    snprintf(dset_name, sizeof(dset_name), "%s%s", qname, "4d_im");
    ReadDataParallelHDF5(rank, dimsf, count, offset, 
            fname, gname, dset_name, data);
    for (i=0; i<len; i++) {
        cin[i][1] = data[i];
    }

    /* Execute and destroy the FFTW plan. */
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    /* Write velocity to HDF5 file */
    int newfile_flag=0, newgroup_flag=0;
    if (strstr(qname, "ux")) {
        newfile_flag = 1;
        newgroup_flag = 1;
    }
    fname = "u4d.h5";
    gname = "/u4d";
    snprintf(dset_name, sizeof(dset_name), "%s%s", qname, "4d");
    dimsf[0] = nt; dimsf[1] = nz;
    dimsf[2] = ny; dimsf[3] = (nx/2+1)*2;
    count[0] = local_n0; count[1] = nz;
    count[2] = ny; count[3] = (nx/2+1)*2;
    offset[0] = local_0_start;
    offset[1] = 0; offset[2] = 0; offset[3] = 0;
    WriteDataParallelHDF5(rank, dimsf, count, offset, fname, 
            gname, dset_name, newfile_flag, newgroup_flag, 1, rout);

    /* Get the maximum value of the data set. */
    double dmax = 0.0; /* Maximum of the data */
    len *= 2; // local_n0*nz*ny*(nx/2+1)*2
    for (i=0; i<len; i++) {
        if (rout[i]>dmax) dmax = rout[i];
    }
    double *dmax_array = (double *)malloc(sizeof(double)*num_procs);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&dmax, 1, MPI_DOUBLE, dmax_array,
            1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    *data_max = 0.0;
    if (my_id == 0) {
        for (i=0; i<num_procs; i++) {
            if (abs(dmax_array[i]>*data_max)) {
                *data_max = abs(dmax_array[i]);
            }
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(data_max, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    fftw_free(data);
    fftw_free(cin);
    fftw_free(rout);
    free(dmax_array);
}

/******************************************************************************
 * Normalize the data to the maximum of the data sets.
 *
 * Input:
 *  nx, ny, nz, nt: the dimensions of the data sets.
 *  vmax: maximum of the all of the 3 velocity components.
 *  qname: quantity name, ux, uy, uz
 ******************************************************************************/
void normalization(int nx, int ny, int nz, int nt, double vmax, char *qname)
{
    int rank = 4;
    const ptrdiff_t n[4] = {nt, nz, ny, (nx/2+1)*2};
    ptrdiff_t alloc_local, local_n0, local_0_start;
    char *fname, *gname, dset_name[256];
    double *data;
    ptrdiff_t i, len;

    /* Get local data size and allocate.
     * This fucntion is from FFTW package.*/
    alloc_local = fftw_mpi_local_size(4, n, MPI_COMM_WORLD, 
            &local_n0, &local_0_start);
    fname = "u4d.h5";
    gname = "/u4d";
    snprintf(dset_name, sizeof(dset_name), "%s%s", qname, "4d");
    hsize_t dimsf[rank], count[rank], offset[rank];
    dimsf[0] = nt; dimsf[1] = nz;
    dimsf[2] = ny; dimsf[3] = (nx/2+1)*2;
    count[0] = local_n0; count[1] = nz;
    count[2] = ny; count[3] = (nx/2+1)*2;
    offset[0] = local_0_start;
    offset[1] = 0; offset[2] = 0; offset[3] = 0;
    data = fftw_alloc_real(alloc_local);
    ReadDataParallelHDF5(rank, dimsf, count, offset, 
            fname, gname, dset_name, data);
    len = local_n0*nz*ny*(nx/2+1)*2;
    for (i=0; i<len; i++) {
        data[i] /= vmax;
    }
    long int sz = len;
    char buff[256];
    snprintf(buff, sizeof(buff), "%s%s", qname, "_distribution");
    GetLinearDistribution(sz, -1.0, 1.0, 100, buff, data);
    GetKineticEnergyDensity(local_n0, nt, nz, ny, (nx/2+1)*2, data);
    WriteDataParallelHDF5(rank, dimsf, count, offset, fname, 
            gname, dset_name, 0, 0, 0, data);
    fftw_free(data);
}

/******************************************************************************
 * Get the distribution of the data in an 1D array using linear intervals.
 *
 * Input:
 *  sz: the size of the array.
 *  dmin, dmax: the minimum and the maximum of the data.
 *  nbands: number of bins for the distribution.
 *  qname: the quantity name.
 *  data: the 1D data array.
 ******************************************************************************/
void GetLinearDistribution(long int sz, double dmin, double dmax,
        int nbands, char *qname, double *data)
{
    double interval;
    long int i;
    int ibin, j;
    int *distribution, *distribution_tot;
    double *base_array;
    distribution = (int *)malloc(sizeof(int)*(nbands+1));
    distribution_tot = (int *)malloc(sizeof(int)*(nbands+1));
    base_array = (double *)malloc(sizeof(double)*(nbands+1));
    interval = (dmax-dmin) / nbands;
    for (j=0; j<nbands+1; j++) {
        distribution[j] = 0;
        distribution_tot[j] = 0;
        base_array[j] = dmin + j*interval;
    }
    for (i=0; i<sz; i++) {
        ibin = floor((data[i]-dmin)/interval);
        if (ibin < 0) ibin = 0;
        if (ibin > nbands) ibin = nbands;
        distribution[ibin]++;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Reduce(distribution, distribution_tot, nbands+1, MPI_INT, 
            MPI_SUM, 0, MPI_COMM_WORLD);
    if (my_id==0) {
        FILE *fp;
        char fname[256];
        snprintf(fname, sizeof(fname), "%s%s", qname, ".dat");
        fp = fopen(fname, "w");
        for (j=0; j<nbands+1; j++) {
            fprintf(fp, "%14.6e %d\n", base_array[j], distribution[j]);
        }
        fclose(fp);
    }
    free(distribution);
    free(distribution_tot);
    free(base_array);
}

/******************************************************************************
 * Get the average of kinetic energy density.
 *
 * Input:
 *  data: the 1D data array, which is actually saved the data for a 4D array
 *  in aligned memory.
 *  local_n0, nz, ny, nx: the actually dimensions of the data in 4D.
 *  nt: the total size of the 1st dimension, which is time here.
 ******************************************************************************/
void GetKineticEnergyDensity(int local_n0, int nt, int nz,
        int ny, int nx, double *data)
{
    int it, i;
    long int sz;
    double *kene; // Local kinetic energy array of current MPI process.
    double *kene_all; // Global kinetic energy array.
    double ene_tot;
    sz = nx * ny * nz;
    kene = (double *)malloc(sizeof(double)*local_n0);
    kene_all = (double *)malloc(sizeof(double)*nt);
    for (it=0; it<local_n0; it++) {
        ene_tot = 0.0;
        for (i=0; i<sz; i++) {
            ene_tot += pow(data[it*sz+i], 2);
        }
        kene[it] = ene_tot / sz;
    }
    for (i=0; i<nt; i++) {
        kene_all[i] = 0.0;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(kene, local_n0, MPI_DOUBLE, kene_all,
            local_n0, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    if (my_id==0) {
        FILE *fp;
        fp = fopen("kinetic_energy.dat", "w");
        for (i=0; i<nt; i++) {
            fprintf(fp, "%d %14.6e\n", i, kene_all[i]);
        }
        fclose(fp);
    }
    free(kene);
    free(kene_all);
}
