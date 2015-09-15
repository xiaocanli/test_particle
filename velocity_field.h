/* Velocity fields */
#ifndef VFIELDS_DOUBLE_H
#define VFIELDS_DOUBLE_H
typedef struct vfields_double {
	double vx, vy, vz;
} vfields_double;
#endif

#ifndef VFIELDS_FLOAT_H
#define VFIELDS_FLOAT_H
typedef struct vfields_float {
	float vx, vy, vz;
} vfields_float;
#endif

void set_variables_velocity(grids *sgrid, domain *sdomain, double v0_field,
        int mt, int dtype);

void initialize_vfield(void);

void free_vfield(void);

void read_vfields_h5(int ct, char *fname, char *gname, int data_type);

void read_vfields_binary(char *filepath, int ct, int data_type);

void get_double_velocity_at_point(double x, double y, double z, double t,
        double *vx, double *vy, double *vz);

void get_float_velocity_at_point(double x, double y, double z, double t,
        float *vx, float *vy, float *vz);
