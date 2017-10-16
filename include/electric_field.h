/* Electric field */
#ifndef EFIELDS_DOUBLE_H
#define EFIELDS_DOUBLE_H
typedef struct efields_double {
	double ex, ey, ez;
} efields_double;
#endif

#ifndef EFIELDS_FLOAT_H
#define EFIELDS_FLOAT_H
typedef struct efields_float {
	float ex, ey, ez;
} efields_float;
#endif

void set_variables_efield(grids *sgrid, domain *sdomain, double B0_field,
        int mt, int dtype);

void initialize_efield(void);

void free_efield(void);

void read_efields_binary(char *filepath, int ct, int data_type);

void get_double_efield_at_point(double x, double y, double z, double t,
        double *ex, double *ey, double *ez);

void get_float_efield_at_point(double x, double y, double z, double t,
        float *ex, float *ey, float *ez);
