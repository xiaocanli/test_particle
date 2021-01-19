#ifndef EMFIELDS_H
#define EMFIELDS_H

#include "hdf5.h"

/* Electromagnetic fields */
typedef struct emfields {
	double Bx, By, Bz;
	double Ex, Ey, Ez;
} emfields;

#ifndef BFIELDS_H
#define BFIELDS_H
/* Magnetic fields */
typedef struct bfields {
	double Bx, By, Bz;
} bfields;
#endif

#ifndef EFIELDS_H
#define EFIELDS_H
/* Electric fields */
typedef struct efields {
	double Ex, Ey, Ez;
} efields;
#endif

void get_emf(double x, double y, double z, double t, struct emfields *emf);
void get_magnetic_field(double x, double y, double z, double t, struct emfields *emf);
void create_fields_ctype(hid_t *memtype, hid_t *filetype);
void getemf_test(double x, double y, double z, double t, struct emfields *emf);
void set_variables_emfields(int stype);
void getemf_mhd_test_particle(double x, double y, double z, double t,
        struct emfields *emf);
void getemf_pic_test_particle(double x, double y, double z, double t,
        struct emfields *emf);
void get_bfield_mhd_test_particle(double x, double y, double z, double t,
        struct emfields *emf);
void get_bfield_pic_test_particle(double x, double y, double z, double t,
        struct emfields *emf);

#endif
