#define I0 4.0E8
//#define I0 4.0E7
#define L0 6.96E8
#define c0 3.0E10

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

/* Electromagnetic fields */
typedef struct emfields {
	double Bx, By, Bz;
	double Ex, Ey, Ez;
} emfields;

/* Magnetic fields */
typedef struct bfields {
	double Bx, By, Bz;
} bfields;

/* Electric fields */
typedef struct efields {
	double Ex, Ey, Ez;
} efields;

/* Velocity fields */
typedef struct vfields {
	double vx, vy, vz;
} vfields;

/* Nested structure for particle information output. */
typedef struct ptl_fields {
    struct particles ptl;
    struct emfields emf;
} ptl_fields;

/* Wire-loop current systems */
typedef struct wlcs{
	double x_wr, y_wr, z_wr;
	double costheta_wr, cosphi_wr;
	double sintheta_wr, sinphi_wr;
    double cur_wr, t0;
	double x_lp, y_lp, z_lp, r_lp;
	double cosalpha_lp, cosbeta_lp;
	double sinalpha_lp, sinbeta_lp;
    double cur_lp;
    double omega;
} wlcs;

extern int ierr, num_procs, my_id;
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
