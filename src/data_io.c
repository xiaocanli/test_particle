#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include <hdf5.h>

/******************************************************************************
 * Read data from HDF5 file use one process. The data is in double.
 *
 * Input:
 *  count: the local dimensions of the data in current MPI process.
 *  offset: the local offset in each dimension in current MPI process.
 *  fname: the name of the HDF5 file.
 *  gname: the name of the group.
 *  dset_name: dataset name.
 *
 * Output:
 *  data: the read data from the file.
 ******************************************************************************/
void read_data_serial_double_h5(int rank, hsize_t *count, hsize_t *offset,
        char *fname, char *gname, char *dset_name, double *data)
{
    hid_t file_id, group_id, dset_id;
    hid_t filespace, memspace;

    /* Open the existing file. */
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Open a group */
    group_id = H5Gopen1(file_id, gname);

    /* Open a dataset and get its dataspace */
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    filespace = H5Dget_space(dset_id);

    memspace = H5Screate_simple(rank, count, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dread(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, H5P_DEFAULT, data);

    /* Close/release resources */
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Read data from HDF5 file use one process. The data is in float.
 ******************************************************************************/
void read_data_serial_float_h5(int rank, hsize_t *count, hsize_t *offset,
        char *fname, char *gname, char *dset_name, float *data)
{
    hid_t file_id, group_id, dset_id;
    hid_t filespace, memspace;

    /* Open the existing file. */
    file_id = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
    /* Open a group */
    group_id = H5Gopen1(file_id, gname);

    /* Open a dataset and get its dataspace */
    dset_id = H5Dopen(group_id, dset_name, H5P_DEFAULT);
    filespace = H5Dget_space(dset_id);

    memspace = H5Screate_simple(rank, count, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, count, NULL);
    H5Dread(dset_id, H5T_NATIVE_FLOAT, memspace, filespace, H5P_DEFAULT, data);

    /* Close/release resources */
    H5Dclose(dset_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Read data from a binary file using one MPI process. The data is in float.
 *
 * Input:
 *  fname: the file name.
 *  offset: the offset from the beginning of the file.
 *  size: the size of the data in number of data.
 *
 * Output:
 *  data: the read data from the file.
 ******************************************************************************/
void read_data_serial_float_binary(char *fname, long int offset, long int size,
        float *data)
{
    FILE *fp;
    fp = fopen(fname, "r");
    if (!fp) {
        printf("Unable to open file %s!\n", fname);
        exit(1);
    }
    fseek(fp, offset * sizeof(float), SEEK_SET);
    fread(data, sizeof(float), size, fp);
    fclose(fp);
}

/******************************************************************************
 * Read data from a binary file using one MPI process. The data is in double.
 ******************************************************************************/
void read_data_serial_double_binary(char *fname, long int offset, long int size,
        double *data)
{
    FILE *fp;
    fp = fopen(fname, "r");
    if (!fp) {
        printf("Unable to open file %s!\n", fname);
        exit(1);
    }
    fseek(fp, offset * sizeof(double), SEEK_SET);
    fread(data, sizeof(double), size, fp);
    fclose(fp);
}
