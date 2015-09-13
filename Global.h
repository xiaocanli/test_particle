#include "emfields.h"

#ifndef INF
#define INF 2e10f
#endif

#ifndef LEN_MAX
#define LEN_MAX 80
#endif

#define charge_mass_proton 9.58E3

#define MAX(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
#define MIN(a,b) \
    ({ __typeof__ (a) _a = (a); \
     __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })

#define nbins 102        /* Total energy bins. */
#define nvar 6           /* x,y,z,vx,vy,vz */
#define ntraj_out 10000  /* Number of trajectory points output. */

#define NUM_THREADS 4

//#define MAX(x, y) (((x) > (y)) ? (x) : (y))
//#define MIN(x, y) (((x) < (y)) ? (x) : (y))

typedef enum {false, true} bool;

// Particle information
typedef struct particles{
    double x, y, z;
    double vx, vy, vz, t;
} particles;

/* Domain information for grids */
typedef struct grids{
    int nx, ny, nz, nt;
    double dx, dy, dz, dt;
} grids;

/* Parameters for force-free field */
typedef struct param_ff {
    double B0, lambda0;
} param_ff;

/* Velocity fields */
typedef struct vfields {
	double vx, vy, vz;
} vfields;

#ifndef PTL_FIELDS_H
#define PTL_FIELDS_H
/* Nested structure for particle information output. */
typedef struct ptl_fields {
    struct particles ptl;
    struct emfields emf;
} ptl_fields;
#endif

extern int ierr, mpi_size, mpi_rank;
extern int isystem; /* flag for the system to use */
//extern struct particles *ptl;
extern struct wlcs *config;
extern int nwlcs, iBessel;
extern struct grids simul_grid;
extern struct param_ff param;
extern double pmass, pcharge, charge_mass;
extern int charge_sign;
extern int nt_out;
extern int nptl_tot; /* Total number of particles. */
extern int ntest_ptl_tot; /* Number of test particles */
/* Accumulated # of particles up to current MPI processes. */
extern int *nptl_accumulate;
extern struct vfields *vfd_b, *vfd_a; /* vfield ahead and behind in t */
extern int imultiple; /* flag for whether multiple time slices
                       * are used. 1: multiple, 0: simple.*/
extern int bc_flag; /* Flag for boundary conditions for particles. */
extern int iadaptive; /* Flag for tracking method to use. 
                       * 0 for fixed step, 1 for adaptive step.*/
extern int rerun_flag; /* flag to check if it is a re-run */
extern int *nsteps_ptl_tracking; /* # of steps each particle is tracked. */
extern int nsteps_output; /* The interval for trajectory output. */
extern int *ntraj_accum_global; /* Accumulated trajectory points 
                                 * over MPI processes */
extern int *ntraj_accum; /* Accumulated trajectory points
                                 * over particles in one MPI process. */
extern struct particles *ptl_time; /* Info for time evolution of particles. */

extern void cross_product(double Ax, double Ay, double Az, 
        double Bx, double By, double Bz, 
        double* Cx, double* Cy, double* Cz);
