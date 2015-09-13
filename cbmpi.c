/******************************************************************************
* This file is part of CHAOTICB.
* Copyright (C) <2012-2014> <Xiaocan Li> <xl0009@uah.edu>
* 
* CHAOTICB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.
* 
* CHAOTICB is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
* 
* You should have received a copy of the GNU General Public License
* along with CHAOTICB.  If not, see <http://www.gnu.org/licenses/>.
******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <omp.h>
#include <mpi.h>
/* #include "Global.h" */
#include "constants.h"
#include "domain.h"
/* #include "cbmpi.h" */
/* #include "tracking.h" */
/* #include "diagnostics.h" */
/* #include "force_free.h" */
/* #include "wlcs.h" */

/* struct vfields *vfd_a, *vfd_b; */

/* double pmass, pcharge, charge_mass; */
/* int ierr, mpi_size, mpi_rank; */
/* int isystem; /1* flag for wire-loop current systems *1/ */
/* int charge_sign; */
/* int bc_flag; /1* Flag for boundary conditions for particles. *1/ */
/* int nptl_tot; */
/* int *nptl_accumulate; */
/* int ntest_ptl_tot = 400; */
/* int *nsteps_ptl_tracking; */

/******************************************************************************
 * Main program for test particle simulation.
 ******************************************************************************/
int main(int argc, char **argv)
{
    int ierr, mpi_size, mpi_rank;
    int nptl, ipvd;
    double vthe, dt;
    char config_file_name[LEN_MAX];

    ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ipvd);
    /* ierr = MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &ipvd); */
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    snprintf(config_file_name, "%s", "init.dat");
    struct domain simul_domain;
    read_domain(mpi_rank, config_file_name, &simul_domain);

    /* assign_ptl(&nptl, &vthe); */
    /* struct particles *ptl; */
    /* ptl = (struct particles*)malloc(sizeof(struct particles)*nptl); */
    /* read_rerun_flag(); */
    /* if (rerun_flag == 0) { */
    /*     /1* new run *1/ */
    /*     init_ptl(nptl, vthe, simul_domain, ptl); */
    /* } else { */
    /*     /1* load particles from last run *1/ */
    /*     read_particles(nptl, ptl, "particles_fields_init.h5"); */
    /* } */

    /* /1* Number of steps each particle is tracked. *1/ */
    /* nsteps_ptl_tracking = (int *)malloc(sizeof(int)*nptl); */

    /* dt = 1.0E-5; */
    /* energy_spectrum(nptl, 0, ptl); */
    /* save_particles_fields(nptl, ptl, nptl_tot, */ 
    /*         nptl_accumulate, "particles_fields_init.h5"); */
    /* particle_tracking_hybrid(nptl, dt, 0, ptl); */
    /* trajectory_diagnostics(nptl, dt, ptl); */

    /* release_memory(ptl); */

    ierr = MPI_Finalize();
	return 0;
}

/******************************************************************************
 * Read initial information of particles.
 *
 * Output:
 *  nptl_tot: total number of particles.
 *  vthe: thermal speed of particles.
 ******************************************************************************/
/* void read_ptlinfo(int *nptl_tot, double *vthe) */
/* { */
    /* //sqrt(kT0/m_p), T0 = 1E6 K, CGS unit */
    /* const double vthe_norm = 9.08E6; */
/*     float np, temp; */
/*     char buff[LEN_MAX]; */
/*     FILE *fp; */
/*     int msg; */
/*     fp = fopen("init.dat", "r"); */
/*     while (fgets(buff,LEN_MAX,fp) != NULL){ */
/*         if (strstr(buff, "Particle information")) { */
/*             break; */
/*         } */
/*     }; */
/*     msg = fscanf(fp, "Total number of particles:%e\n", &np); */
/*     if (msg != 1) { */
/*         printf("Failed to read the total number of particles.\n"); */
/*         exit(1); */
/*     } */
/*     *nptl_tot = (int)np; */
/*     msg = fscanf(fp, "Particle temperature:%f\n", &temp); */
/*     if (msg != 1) { */
/*         printf("Failed to read the particle temperature.\n"); */
/*         exit(1); */
/*     } */
/*     msg = fscanf(fp, "Mass:%lf\n", &pmass); */
/*     if (msg != 1) { */
/*         printf("Failed to read the particle mass.\n"); */
/*         exit(1); */
/*     } */
/*     msg = fscanf(fp, "Charge:%lf\n", &pcharge); */
/*     if (msg != 1) { */
/*         printf("Failed to read the particle charge.\n"); */
/*         exit(1); */
/*     } */
/*     charge_sign = pcharge/fabs(pcharge); */
/*     *vthe = sqrt(temp/pmass)*vthe_norm; */
/*     charge_mass = charge_mass_proton*pcharge/pmass; */

    // Boundary condition for particles
/*     while (fgets(buff,LEN_MAX,fp) != NULL){ */
/*         if (strstr(buff, "0: periodic 1: open")) { */
/*             break; */
/*         } */
/*     }; */
/*     msg = fscanf(fp, "BC for particles: %d\n", &bc_flag); */
/*     if (msg != 1) { */
/*         printf("Failed to read BC for particles.\n"); */
/*         exit(1); */
/*     } */

/*     if (mpi_rank == 0) { */
/*         printf("========================================\n"); */
/*         printf("Total number of particles = %e\n", np); */
/*         printf("Initial temperature = %f\n", temp); */
/*         printf("Mass of particles = %f\n", pmass); */
/*         printf("Charge of particles = %f\n", pcharge); */
/*         printf("Particle thermal velocity = %f\n", *vthe); */
/*         if (bc_flag == 0) { */
/*             printf("Using periodic boundary condition for particles\n"); */
/*         } */
/*         else if (bc_flag == 1) { */
/*             printf("Using open boundary condition for particles\n"); */
/*         } */
/*         else { */
/*             printf("ERROR: Wrong boundary conditions."); */
/*         } */
/*         printf("========================================\n"); */
/*     } */
/*     fclose(fp); */
/* } */

/******************************************************************************
 * Generating 2 independent Gaussian random numbers using Box-Muller method
 * Source: section 7.3.4 of Numerical Recipes 2007
 *
 * Output:
 *  ran1, ran2: two Gaussian random numbers.
 ******************************************************************************/
/* void gaussian_rand(double *ran1, double *ran2) */
/* { */
/*     double v1, v2, rsq, ran, fac; */
/*     rsq = 2.0; */
/* //    srand48(time(NULL)); */
/*     while ((rsq >= 1.0) | (rsq == 0.0)) { */
/*         ran = drand48(); */
/*         //printf("%f \n", ran); */
/*         v1 = 2.0*ran - 1.0; */
/*         ran = drand48(); */
/*         //printf("%f \n", ran); */
/*         v2 = 2.0*ran - 1.0; */
/*         rsq = v1*v1 + v2*v2; */
/*     } */
/*     fac = sqrt(-2.0*log(rsq)/rsq); */
/*     *ran1 = v1*fac; */
/*     *ran2 = v2*fac; */
/* } */

/******************************************************************************
 * Initialization of particles for different systems.
 *
 * Input:
 *  nptl: particle number for this MPI process.
 *  vthe: thermal speed of particles.
 *  simul_domain: simulation domain information.
 * Output:
 *  ptl: particle structure array.
 ******************************************************************************/
/* void init_ptl(int nptl, double vthe, domain simul_domain, struct particles *ptl) */
/* { */
/*     double xmin, xmax, ymin, ymax, zmin, zmax; */
/*     double ran1, ran2, tmp; */
/*     int i; */
/*     xmin = simul_domain.xmin; */
/*     xmax = simul_domain.xmax; */
/*     ymin = simul_domain.ymin; */
/*     ymax = simul_domain.ymax; */
/*     zmin = simul_domain.zmin; */
/*     zmax = simul_domain.zmax; */
/*     tmp = sqrt(2.0) * vthe; */
/*     //printf("Thermal speed: %f\n", vthe); */
/*     srand48((unsigned)time(NULL)+mpi_rank*mpi_size); */
/*     if (isystem == 0) { */
/*         double nvthe = 100.0; */
/*         for (i = 0; i < nptl; i++) { */
/*             //ptl[i].vx = drand48()*nvthe*tmp / c0; */
/*             //ptl[i].vy = drand48()*nvthe*tmp / c0; */
/*             ptl[i].x = 0.0; */
/*             ptl[i].y = 0.0; */
/*             ptl[i].z = 0.0; */
/*             ptl[i].t = 0.0; */
/*             ptl[i].vx = nvthe*(i+1)*tmp / (c0*nptl); */
/*             ptl[i].vy = nvthe*(i+1)*tmp / (c0*nptl); */
/*             ptl[i].vz = 0.0; */
/*             //printf("%f %f %f\n", ptl[i].vx, ptl[i].vy, ptl[i].vz); */
/*         } */
/*     } */ 
/*     else { */
/*         for (i = 0; i < nptl; i++) { */
/*             ptl[i].x = (xmax-xmin)*drand48() + xmin; */
/*             ptl[i].y = (ymax-ymin)*drand48() + ymin; */
/*             ptl[i].z = (zmax-zmin)*drand48() + zmin; */
/*             //printf("%f %f %f\n", ptl[i].x, ptl[i].y, ptl[i].z); */
/*             ptl[i].t = 0.0; */
/*             gaussian_rand(&ran1, &ran2); */
/*             ptl[i].vx = ran1*tmp / c0; */
/*             ptl[i].vy = ran2*tmp / c0; */
/*             //printf("%f %f \n", ran1, ran2); */
/*             gaussian_rand(&ran1, &ran2); */
/*             ptl[i].vz = ran1*tmp / c0; */
/*             //printf("%f %f %f\n", ptl[i].vx, ptl[i].vy, ptl[i].vz); */
/*         } */
/*     } */
/* } */

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
 * Read particle initial information, assign particles to each process and
 * initialize the particles.
 *
 * Output:
 *  nptl: total number of particles of each MPI process.
 *  vthe: thermal speed of particles.
 ******************************************************************************/
/* void assign_ptl(int *nptl, double *vthe) */
/* { */
/*     /1* Read from file the total particle number and particle thermal speed *1/ */
/*     read_ptlinfo(&nptl_tot, vthe); */
/*     if (mpi_rank == 0) { */
/*         printf("Tracking %d particles using %d processes\n", */ 
/*                 nptl_tot, mpi_size); */
/*     } */

/*     /1* Number of test particles for trajectory diagnostics. *1/ */
/*     if (ntest_ptl_tot > nptl_tot) { */
/*         ntest_ptl_tot = nptl_tot; */
/*     } */

/*     if (isystem == 0) { */
/*         nptl_tot = ntest_ptl_tot; */
/*     } */
/*     nptl_accumulate = (int *)malloc(sizeof(int)*mpi_size); */
/*     particle_broadcast(nptl_tot, nptl); */
/*     printf("MPI process %d will trace %d particles. \n", mpi_rank, *nptl); */
/* } */
