/* Velocity fields */
#ifndef BFIELDS_DOUBLE_H
#define BFIELDS_DOUBLE_H
typedef struct bfields_double {
	double bx, by, bz;
} bfields_double;
#endif

#ifndef BFIELDS_FLOAT_H
#define BFIELDS_FLOAT_H
typedef struct bfields_float {
	float bx, by, bz;
} bfields_float;
#endif

void set_variables_bfield(grids *sgrid, domain *sdomain, double B0_field,
        int mt, int dtype);

void initialize_bfield(void);

void free_bfield(void);

void read_bfields_h5(int ct, char *fname, char *gname, int data_type);

void read_bfields_binary(char *filepath, int ct, int data_type);

void get_double_bfield_at_point(double x, double y, double z, double t,
        double *bx, double *by, double *bz);

void get_float_bfield_at_point(double x, double y, double z, double t,
        float *bx, float *by, float *bz);
