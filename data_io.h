#include "hdf5.h"

void read_data_serial_double_h5(int rank, hsize_t *count, hsize_t *offset,
        char *fname, char *gname, char *dset_name, double *data);

void read_data_serial_float_h5(int rank, hsize_t *count, hsize_t *offset,
        char *fname, char *gname, char *dset_name, float *data);

void read_data_serial_double_binary(char *fname, long int offset, long int size,
        double *data);

void read_data_serial_float_binary(char *fname, long int offset, long int size,
        float *data);
