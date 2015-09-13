#ifndef EMFIELDS_H
#define EMFIELDS_H
/* Electromagnetic fields */
typedef struct emfields {
	double Bx, By, Bz;
	double Ex, Ey, Ez;
} emfields;
#endif

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
