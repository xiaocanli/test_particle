#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>
#include "hdf5.h"
#include "constants.h"
#include "particle_info.h"
#include "domain.h"

int read_rerun_flag(char *config_file_name);
void init_ptl(int mpi_rank, int mpi_size, domain simul_domain, int nptl,
        double vthe, int system_type, struct particles *ptl);

/******************************************************************************
 * Assign particles to each process.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  mpi_size: total number of MPI processes.
 *  system_type: the type of system.
 *
 * Output:
 *  nptl: total number of particles of each MPI process.
 *
 * Input & output:
 *  nptl_tot: total number of particles.
 *  nptl_traj_tot: number of particles for trajectory diagnostics.
 ******************************************************************************/
void assign_particles(int mpi_rank, int mpi_size, int *nptl_tot,
        int *nptl_traj_tot, int system_type, int *nptl, int *nptl_accumulate)
{
    /* Number of test particles for trajectory diagnostics. */
    if ((*nptl_traj_tot) > (*nptl_tot)) {
        *nptl_traj_tot = (*nptl_tot);
    }

    if (system_type == 0) {
        /* test system. */
        *nptl_tot = (*nptl_traj_tot);
    }

    particle_broadcast(mpi_rank, mpi_size, *nptl_tot, nptl, nptl_accumulate);
    if (mpi_rank == 0) {
        printf("MPI process %d will trace %d particles. \n",
                mpi_rank, *nptl);
    }
}

/******************************************************************************
 * Assign particles to MPI process and get the array for accumulated particle
 * number from MPI process 0 to mpi_size.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  mpi_size: total number of MPI processes.
 *  nptl_tot: total # of particles.
 *
 * Output:
 *  nptl: # of particles current MPI process.
 *  nptl_accumulate: global array for accumulated particle number.
 ******************************************************************************/
void particle_broadcast(int mpi_rank, int mpi_size, int nptl_tot, int *nptl,
        int *nptl_accumulate)
{
    int i;
    if (mpi_rank == 0) {
        for (i = 0; i < mpi_size; i++) {
            /* Number of particles for each process. */
            nptl_accumulate[i] = nptl_tot / mpi_size;
            if (i < (nptl_tot % mpi_size)) {
                nptl_accumulate[i] += 1;
            }
        }
        for (i = 1; i < mpi_size; i++) {
            /* Accumulation to get the shift. */
            nptl_accumulate[i] += nptl_accumulate[i-1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(nptl_accumulate, mpi_size, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        *nptl = nptl_accumulate[0];
    }
    else {
        *nptl = nptl_accumulate[mpi_rank]-nptl_accumulate[mpi_rank-1];
    }
}

/******************************************************************************
 * Read initial information of particles.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  config_file_name: the configuration filename.
 *
 * Output:
 *  nptl_tot: total number of particles.
 *  ptemp: particle temperature.
 *  pmass: particle mass.
 *  pcharge: particle charge.
 *  charge_sign: the sign of the particles.
 *  vthe: thermal speed of particles.
 *  charge_mass: the ratio of particle charge and mass.
 *  bc_flag: the boundary condition flag.
 *  nptl_traj_tot: number of particles for trajectory diagnostics.
 ******************************************************************************/
void read_particle_info(int mpi_rank, char *config_file_name, int *nptl_tot,
        double *ptemp, double *pmass, double *pcharge, int *charge_sign,
        double *vthe, double *charge_mass, int *bc_flag, int *nptl_traj_tot)
{
    float np;
    char buff[LEN_MAX];
    FILE *fp;
    int msg;
    fp = fopen(config_file_name, "r");
    while (fgets(buff,LEN_MAX,fp) != NULL){
        if (strstr(buff, "Particle information")) {
            break;
        }
    };
    msg = fscanf(fp, "Total number of particles:%e\n", &np);
    if (msg != 1) {
        printf("Failed to read the total number of particles.\n");
        exit(1);
    }
    *nptl_tot = (int)np;
    msg = fscanf(fp, "Particle temperature:%lf\n", ptemp);
    if (msg != 1) {
        printf("Failed to read the particle temperature.\n");
        exit(1);
    }
    msg = fscanf(fp, "Mass:%lf\n", pmass);
    if (msg != 1) {
        printf("Failed to read the particle mass.\n");
        exit(1);
    }
    msg = fscanf(fp, "Charge:%lf\n", pcharge);
    if (msg != 1) {
        printf("Failed to read the particle charge.\n");
        exit(1);
    }

    msg = fscanf(fp, "Number of particles for trajectory diagnostics:%d\n",
            nptl_traj_tot);
    if (msg != 1) {
        printf("Failed to read the number of particles for trajectory diagnostics\n");
        exit(1);
    }

    *charge_sign = (*pcharge) / fabs(*pcharge);
    *vthe = sqrt((*ptemp) / (*pmass)) * VTHE_NORM;
    *charge_mass = CHARGE_MASS_PROTON * (*pcharge) / (*pmass);

    /* Boundary condition for particles*/
    while (fgets(buff, LEN_MAX, fp) != NULL){
        if (strstr(buff, "0: periodic 1: open")) {
            break;
        }
    };
    msg = fscanf(fp, "BC for particles: %d\n", bc_flag);
    if (msg != 1) {
        printf("Failed to read BC for particles.\n");
        exit(1);
    }

    if (mpi_rank == 0) {
        printf("========================================\n");
        printf("Total number of particles = %e\n", (double)*nptl_tot);
        printf("Initial temperature = %f\n", *ptemp);
        printf("Mass of particles = %f\n", *pmass);
        printf("Charge of particles = %f\n", *pcharge);
        printf("Particle thermal velocity = %f\n", *vthe);
        if ((*bc_flag) == 0) {
            printf("Using periodic boundary condition for particles\n");
        }
        else if ((*bc_flag) == 1) {
            printf("Using open boundary condition for particles\n");
        }
        else {
            printf("ERROR: Wrong boundary conditions.");
        }
        printf("========================================\n");
    }
    fclose(fp);
}

/******************************************************************************
 * Initialize particles either analytically or loading from previous runs.
 * 
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  mpi_size: total number of MPI processes.
 *  config_file_name: the configuration filename.
 *  nptl: number of particles for current MPI process.
 *  vthe: thermal speed of particles.
 *  system_type: the system type.
 *
 * Output:
 *  ptl: particle structure array.
 ******************************************************************************/
void initialize_partilces(int mpi_rank, int mpi_size, char *config_file_name,
        domain simul_domain, int nptl, double vthe, int system_type,
        int *nptl_accumulate,  particles *ptl)
{
    int rerun_flag;
    rerun_flag = read_rerun_flag(config_file_name);
    if (mpi_rank == 0) {
        if (rerun_flag == 1) {
            printf("This run is going to rerun the last simulation.\n");
        } else {
            printf("This is a new simulation.\n");
        }
    }
    if (rerun_flag == 0) {
        /* new run */
        init_ptl(mpi_rank, mpi_size, simul_domain, nptl, vthe, system_type, ptl);
    } else {
        /* load particles from last run */
        read_particles(mpi_rank, nptl, nptl_accumulate, ptl,
                "particles_fields_init.h5");
    }
}

/******************************************************************************
 * Read the flag to check if it is a re-run of previous simulation.
 * 
 * Return:
 *  rerun_flag: 0 for new run. 1 for re-run.
 ******************************************************************************/
int read_rerun_flag(char *config_file_name)
{
    FILE *fp;
    char buff[LEN_MAX];
    int msg, rerun_flag;
    fp = fopen(config_file_name, "r");
    while (fgets(buff,LEN_MAX,fp) != NULL){
        if (strstr(buff, "Check points info")) {
            break;
        }
    };
    msg = fscanf(fp, "re-run the simulation:%d\n", &rerun_flag);
    if (msg != 1) {
        printf("Failed to read the re-run flag.\n");
        exit(1);
    }
    fclose(fp);
    return rerun_flag;
}

/******************************************************************************
 * Generating 2 independent Gaussian random numbers using Box-Muller method
 * Source: section 7.3.4 of Numerical Recipes 2007
 *
 * Output:
 *  ran1, ran2: two Gaussian random numbers.
 ******************************************************************************/
void gaussian_rand(double *ran1, double *ran2)
{
    double v1, v2, rsq, ran, fac;
    rsq = 2.0;
    while ((rsq >= 1.0) | (rsq == 0.0)) {
        ran = drand48();
        v1 = 2.0*ran - 1.0;
        ran = drand48();
        v2 = 2.0*ran - 1.0;
        rsq = v1*v1 + v2*v2;
    }
    fac = sqrt(-2.0*log(rsq)/rsq);
    *ran1 = v1*fac;
    *ran2 = v2*fac;
}

/******************************************************************************
 * Initialization of particles for different systems.
 *
 * Input:
 *  nptl: particle number for this MPI process.
 *  vthe: thermal speed of particles.
 * Output:
 *  ptl: particle structure array.
 ******************************************************************************/
void init_ptl(int mpi_rank, int mpi_size, domain simul_domain, int nptl,
        double vthe, int system_type, struct particles *ptl)
{
    double xmin, xmax, ymin, ymax, zmin, zmax;
    double ran1, ran2, tmp;
    int i;
    xmin = simul_domain.xmin;
    xmax = simul_domain.xmax;
    ymin = simul_domain.ymin;
    ymax = simul_domain.ymax;
    zmin = simul_domain.zmin;
    zmax = simul_domain.zmax;
    tmp = sqrt(2.0) * vthe;
    //printf("Thermal speed: %f\n", vthe);
    srand48((unsigned)time(NULL)+mpi_rank*mpi_size);
    if (system_type == 0) {
        double nvthe = 100.0;
        for (i = 0; i < nptl; i++) {
            //ptl[i].vx = drand48()*nvthe*tmp / c0;
            //ptl[i].vy = drand48()*nvthe*tmp / c0;
            ptl[i].x = 0.0;
            ptl[i].y = 0.0;
            ptl[i].z = 0.0;
            ptl[i].t = 0.0;
            ptl[i].vx = nvthe*(i+1)*tmp / (c0*nptl);
            ptl[i].vy = nvthe*(i+1)*tmp / (c0*nptl);
            ptl[i].vz = 0.0;
            //printf("%f %f %f\n", ptl[i].vx, ptl[i].vy, ptl[i].vz);
        }
    } 
    else {
        for (i = 0; i < nptl; i++) {
            ptl[i].x = (xmax-xmin)*drand48() + xmin;
            ptl[i].y = (ymax-ymin)*drand48() + ymin;
            ptl[i].z = (zmax-zmin)*drand48() + zmin;
            //printf("%f %f %f\n", ptl[i].x, ptl[i].y, ptl[i].z);
            ptl[i].t = 0.0;
            gaussian_rand(&ran1, &ran2);
            ptl[i].vx = ran1*tmp / c0;
            ptl[i].vy = ran2*tmp / c0;
            //printf("%f %f \n", ran1, ran2);
            gaussian_rand(&ran1, &ran2);
            ptl[i].vz = ran1*tmp / c0;
            //printf("%f %f %f\n", ptl[i].vx, ptl[i].vy, ptl[i].vz);
        }
    }
}

/******************************************************************************
 * Calculate electromagnetic fields for test case.
 * 
 * Input:
 * x, y, z, t: spatial positions and time.
 * Output:
 * emf: electromagnetic fields at (x,y,z,t) 
 ******************************************************************************/
/* void getemf_test(double x, double y, double z, double t, struct emfields *emf) */
/* { */
/*     emf->Bx = 0.0; */
/*     emf->By = 0.0; */
/*     emf->Bz = 1.0E0; */
/*     emf->Ex = 0.0; */
/*     emf->Ey = 0.0; */
/*     emf->Ez = 0.0; */
/* } */

/******************************************************************************
 * Calculate or interpolate to get electromagnetic fields.
 * isystem = 0: test case.
 * isystem = 1: wire-loop current system case.
 * isystem = 2: force-free magnetic field case.
 *
 * Input:
 *  x, y, z, t: spatial positions and time.
 * Output:
 *  emf: electromagnetic fields at (x,y,z,t) 
 ******************************************************************************/
/* void get_emf(double x, double y, double z, double t, struct emfields *emf) */
/* { */
/*     if (isystem == 0) { */
/*         getemf_test(x, y, z, t, emf); */
/*     } */
/*     else if (isystem == 1) { */
/*         getemf_wlcs(x, y, z, t, emf); */
/*     } */
/*     else if (isystem == 2){ */
/*         getemf_ff(x, y, z, t, emf); */
/*     } */
/*     //printf("Bx, By, Bz, Ex, Ey, Ez: %10.6e %10.6e %10.6e %10.6e %10.6e %10.6e\n", */
/*     //        emf->Bx, emf->By, emf->Bz, emf->Ex, emf->Ey, emf->Ez); */
/* } */

/******************************************************************************
 * Release allocated memory.
 * 
 * Released memory:
 *  ptl: particle structure array.
 *  config: configuration for wlcs system when isystem = 1.
 ******************************************************************************/
/* void release_memory(struct particles *ptl) */
/* { */
/*     free(ptl); */
/*     if (isystem == 1) { */
/*         free(config); */
/*     } */
/*     else if (isystem == 2) { */
/*         free(vfd_b); */
/*         if (imultiple == 1) { */
/*             free(vfd_a); */
/*         } */
/*     } */
/*     free(nptl_accumulate); */
/*     free(nsteps_ptl_tracking); */
/* } */

/******************************************************************************
 * Read all particles information from hdf5 file.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nptl: number of particle in one MPI process.
 *  nptl_accumulate: global array for accumulated particle number.
 *
 * Output:
 *  ptl: particle structure array.
 ******************************************************************************/
void read_particles(int mpi_rank, int nptl, int *nptl_accumulate,
        struct particles *ptl, char *fname)
{
    int rank = 1;
//    char fname[] = "particles_fields.h5";
    char gname[] = "/particles_fields";
    hid_t file_id, group_id;
    hid_t dset_ptl;
    hid_t filespace, memspace;
    hid_t plist_id;
    //herr_t status;
    hsize_t count[rank], offset[rank];
    hid_t memtype_ptl, filetype_ptl;

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

    /* Create a new file collectively and release property list identifier. */
    file_id = H5Fopen(fname, H5F_ACC_RDWR, plist_id);
    H5Pclose(plist_id);
    /* Create a new group */
    group_id = H5Gopen1(file_id, gname);

    /* Create compound datatype */
    create_particles_ctype(&memtype_ptl, &filetype_ptl);

    /* Open a dataset and get its dataspace. */
    dset_ptl = H5Dopen(group_id, "particles", H5P_DEFAULT);
    filespace = H5Dget_space(dset_ptl);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    /* Count and offset in the memory. */
    count[0] = nptl;
    if (mpi_rank == 0) {
        offset[0] = 0;
    } else {
        offset[0] = nptl_accumulate[mpi_rank-1];
    }

    memspace = H5Screate_simple(rank, count, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, 
            NULL, count, NULL);
    H5Dread(dset_ptl, memtype_ptl, memspace, filespace, 
            plist_id, ptl);
    H5Sclose(filespace);

    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Dclose(dset_ptl);
    H5Tclose(memtype_ptl);
    H5Tclose(filetype_ptl);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Create a compound data type containing the particle information for HDF5.
 *
 * Output:
 *  memtype: compound datatype for memory.
 *  filetype: compound datatype for the file.
 ******************************************************************************/
void create_particles_ctype(hid_t *memtype, hid_t *filetype)
{
    //herr_t status;
    /* Create the compound datatype in memory */
    *memtype = H5Tcreate(H5T_COMPOUND, sizeof(particles));
    H5Tinsert(*memtype, "x", 
            HOFFSET(particles, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "y", 
            HOFFSET(particles, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "z", 
            HOFFSET(particles, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "vx", 
            HOFFSET(particles, vx), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "vy", 
            HOFFSET(particles, vy), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "vz", 
            HOFFSET(particles, vz), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "t", 
            HOFFSET(particles, t), H5T_NATIVE_DOUBLE);
    /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */
    *filetype = H5Tcreate(H5T_COMPOUND, 8*7);
    H5Tinsert (*filetype, "x", 0, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "y", 8, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "z", 8*2, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "vx", 8*3, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "vy", 8*4, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "vz", 8*5, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "t", 8*6, H5T_IEEE_F64BE);
}
