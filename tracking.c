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
#include "constants.h"
#include "domain.h"
#include "tracking.h"
#include "diagnostics.h"
#include "StepperBS.h"
#include "particle_info.h"
#include "domain.h"
#include "emfields.h"

struct grids simul_grid;
struct domain simul_domain;
int system_type;
int multi_tframe;
int bc_flag;
double charge_mass;
int charge_sign;

void init_spectrum(int nbins, int nt_out, double espectrum [][nbins]);

/******************************************************************************
 * Set shared variables for current file.
 ******************************************************************************/
void set_variables_tracking(grids *sgrid, domain *sdomain, int stype, int mt)
{
    memcpy(&simul_grid, sgrid, sizeof(grids));
    memcpy(&simul_domain, sdomain, sizeof(domain));
    system_type = stype;
    multi_tframe = mt;
}

void set_ptl_params_tracking(int pbc, double cmass, int csign)
{
    bc_flag = pbc;
    charge_mass = cmass;
    charge_sign = csign;
}

/******************************************************************************
 * Particle tracking using fixed time steps.
 *
 * Input:
 *  nptl: particle number for this process.
 *  dt: time step for particle tracking.
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *  traj_diagnose: flag for whether to do particle trajectory diagnostics.
 *  nsteps_ptl_tracking: the number of tracking steps for each particle.
 *
 * Output:
 *  ptl: particles structure array will be updated.
 *  espectrum, espect_escape: spectra summed over all threads.
 *  espect_private, espect_escape_private: spectra in each thread.
 *  espect_tot: spectra summed over all MPI process.
 ******************************************************************************/
void particle_tracking_fixed(int nptl, double dt, int nbins, int nt_out,
        int traj_diagnose, int nsteps_output, double pmass,
        int *nsteps_ptl_tracking, particles *ptl, int *ntraj_accum,
        double espectrum[][nbins], double espect_tot[][nbins],
        double espect_escape[][nbins], double espect_private[][nbins*nt_out],
        double espect_escape_private[][nbins*nt_out], particles *ptl_time)
{
    int estep, einterval_t, ptl_id;
    int *ntraj_diagnostics_points_array;
    estep = simul_domain.tmax /dt;
    einterval_t = estep/nt_out;
    printf("Energy output time interval %d\n", einterval_t);

    /* Number of steps each particle is tracked. */
    for (ptl_id=0; ptl_id<nptl; ptl_id++) {
        nsteps_ptl_tracking[ptl_id] = 1;
    }

    if (traj_diagnose == 1) {
        /* Temporary array for trajectory diagnostics. 
         * It keeps tracking how many trajectory diagnostics points
         * for each test particle. Doing it as an array because some
         * particles are going to get out the simulation box. */
        ntraj_diagnostics_points_array = (int *)malloc(sizeof(int)*nptl);
        for (ptl_id = 0; ptl_id < nptl; ptl_id++) {
            /* At least one point. */
            ntraj_diagnostics_points_array[ptl_id] = 1;
        }
    }

    #pragma omp parallel
    {
        int j;
        int tid = omp_get_thread_num();
        int num_threads = omp_get_num_threads();
        if (tid == 0) {
            printf("Number of threads: %d\n", num_threads);
        }
        for (j=0; j<nbins*nt_out; j++) {
            espect_private[tid][j] = 0.0;
            espect_escape_private[tid][j] = 0.0;
        }
        if ((system_type == 2 || system_type == 3) && multi_tframe == 1) {
            /* Force-free system and interpolation between multiple
             * velocity field. */
            /* ParticleTrackingOmpFF(nptl, dt, tid, estep, einterval_t, */
            /*         traj_diagnose, ptl, nsteps_ptl_tracking, */
            /*         espect_private, espect_escape_private); */
        } else {
            MultiStepsParticleTrackingOmp(nptl, dt, tid, 1, estep, nbins,
                    nt_out, einterval_t, traj_diagnose, nsteps_output, pmass,
                    ntraj_accum, ntraj_diagnostics_points_array, ptl,
                    nsteps_ptl_tracking, espect_private, 
                    espect_escape_private, ptl_time);
        }
        GatherParticleSpectra(num_threads, nbins, nt_out, espectrum, 
                espect_private, espect_escape, espect_escape_private);
    }

    if (traj_diagnose == 1) {
        free(ntraj_diagnostics_points_array);
    }
}

/******************************************************************************
 * Particle tracking procedure parallelized with OpenMP for force-free field
 * case. The difference here is that we need velocity field to calculate
 * electric field in force-free field case. So we need to update the velocity
 * field to include of the dynamical effects.
 *
 * Input:
 *  nptl: particle number for this process.
 *  dt: time step for particle tracking.
 *  tid: thread ID.
 *  estep: total time steps for particle tracking.
 *  einterval_t: time interval for energy spectra diagnostics.
 *  traj_diagnose: flag for whether to do particle trajectory diagnostics.
 *
 * Input & Output:
 *  ptl: particles struct array will be updated.
 *  nsteps_ptl_tracking: total tracking steps up to now for each particle.
 *  espect_private, espect_escape_private: particle energy spectra inside the
 *      box and escaping the box for current thread.
 ******************************************************************************/
/* void ParticleTrackingOmpFF(int nptl, double dt, int tid, int estep, */
/*         int einterval_t, int traj_diagnose, struct particles *ptl, */
/*         int *nsteps_ptl_tracking, double espect_private[][nbins*nt_out], */ 
/*         double espect_escape_private[][nbins*nt_out]) */
/* { */
/*     double vfield_dt; // The time interval between slices of velocity field. */
/*     int current_slice; // The time slice of velocity field currently in use. */
/*     int starting, ending; */
/*     vfield_dt = simul_domain.tmax / (simul_grid.nt-1); */
/*     current_slice = 0; */
/*     starting = 1; */
/*     ending = floor(((current_slice+1.0)*vfield_dt)/dt); */
/*     for (current_slice = 1; current_slice < simul_grid.nt; current_slice++) { */
/*         MultiStepsParticleTrackingOmp(nptl, dt, tid, starting, ending, */
/*             einterval_t, traj_diagnose, ptl, nsteps_ptl_tracking, */
/*             espect_private, espect_escape_private); */
/*         starting = ending + 1; */
/*         ending = floor(((current_slice+1.0)*vfield_dt)/dt); */
/*         /1* Waiting for all threads to reach this point. *1/ */
/*         #pragma omp barrier */
/*         /1* Only one thread is going to update the velocity field. *1/ */
/*         #pragma omp master */
/*         ReadVelocityField(current_slice); */
/*         /1* Waiting for all threads to reach this point. *1/ */
/*         #pragma omp barrier */
/*     } */
/* } */

/******************************************************************************
 * Read velocity field from file for current time slice.
 * Input:
 *  current_slice: current time slice ID.
 * Output:
 *  vfd_a, vfd_b: global arrays for velocity field behind and ahead of current
 *      time point.
 ******************************************************************************/
/* void ReadVelocityField(int current_slice) */
/* { */
/*     double v0; */
/*     long int sz; */
/*     dim_vfield(&v0); */
/*     sz = simul_grid.nx * simul_grid.ny * simul_grid.nz; */
/*     if (current_slice < simul_grid.nt-1) { */
/*         memcpy(vfd_b, vfd_a, sizeof(struct vfields)*sz); */
/*         read_vfields(current_slice+1, v0, vfd_a); */
/*     } */
/* } */

/******************************************************************************
 * Multiple steps of particle tracking using fixed time step parallelized using
 * OpenMP.
 *
 * Input:
 *  nptl: particle number for this process.
 *  dt: time step for particle tracking.
 *  tid: thread ID.
 *  sstep: starting step for particle tracking.
 *  estep: ending step for particle tracking.
 *  einterval_t: time interval for energy spectra diagnostics.
 *  traj_diagnose: flag for whether to do particle trajectory diagnostics.
 *
 * Input & Output:
 *  ptl: particles struct array will be updated.
 *  nsteps_ptl_tracking: total tracking steps up to now for each particle.
 *  espect_private, espect_escape_private: particle energy spectra inside the
 *      box and escaping the box for current thread.
 ******************************************************************************/
void MultiStepsParticleTrackingOmp(int nptl, double dt, int tid, int sstep,
        int estep, int nbins, int nt_out, int einterval_t, int traj_diagnose,
        int nsteps_output, double pmass, int *ntraj_accum,
        int *ntraj_diagnostics_points_array, struct particles *ptl,
        int *nsteps_ptl_tracking, double espect_private[][nbins*nt_out], 
        double espect_escape_private[][nbins*nt_out], particles *ptl_time)
{
    int ntp, iescape_pre, iescape_after;
    int i, j;
    #pragma omp for private(j,ntp,iescape_pre,iescape_after)
    for (i = 0; i < nptl; i++) {
        iescape_pre = 0; // Default: not escaped yet.
        particle_bc(&ptl[i].x, &ptl[i].y, &ptl[i].z, &iescape_pre);
        iescape_after = iescape_pre;
        ntp = floor(sstep/(double)einterval_t);

        int ntraj_offset_local = 0;
        int ntraj_diagnostics_points = 1;
        if (traj_diagnose == 1) {
            /* Trajectory diagnostics. */
            if (i != 0) {
                ntraj_offset_local = ntraj_accum[i-1];
            }
            ntraj_diagnostics_points = ntraj_diagnostics_points_array[i];
        }

        for (j=sstep; j<=estep; j++) {
            SingleStepParticleTrackingOMP(i, j, dt, tid, einterval_t,
                    traj_diagnose, nbins, nt_out, nsteps_output, pmass, &ntp, ptl,
                    &iescape_pre, &iescape_after, ntraj_offset_local,
                    &ntraj_diagnostics_points, nsteps_ptl_tracking,
                    espect_private, espect_escape_private, ptl_time);
            if (traj_diagnose == 1) {
                ntraj_diagnostics_points_array[i] = ntraj_diagnostics_points;
            }
        }
    }
}

/******************************************************************************
 * A single step of particle tracking using fixed time step parallelized using
 * OpenMP.
 *
 * Input:
 *  iptl: particle index
 *  it: time point of particle tracking.
 *  dt: time step for particle tracking.
 *  ntraj_offset_local: since particle trajectory diagnostics results are
 *      saved in the same array, we need to figure out the starting point
 *      for the data in current thread in this array.
 *  tid: thread ID.
 *  einterval_t: time interval for energy spectra diagnostics.
 *  traj_diagnose: flag for whether to do particle trajectory diagnostics.
 *
 * Input & Output:
 *  ntp: Current # of output time points for energy spectra diagnostics.
 *  ptl: particles struct array will be updated.
 *  iescape_pre, iescape_after: check if particle escape from the box
 *      before and after this step.
 *  ntraj_diagnostics_points: how many time points have been output for particle trajectory
 *      diagnostics.
 *  nsteps_ptl_tracking: total tracking steps up to now for each particle.
 *  espect_private, espect_escape_private: particle energy spectra inside the
 *      box and escaping the box for current thread.
 ******************************************************************************/
void SingleStepParticleTrackingOMP(int iptl, int it, double dt, int tid,
        int einterval_t, int traj_diagnose, int nbins, int nt_out,
        int nsteps_output, double pmass, int *ntp, particles *ptl,
        int *iescape_pre, int *iescape_after, int ntraj_offset_local,
        int *ntraj_diagnostics_points, int *nsteps_ptl_tracking,
        double espect_private[][nbins*nt_out], 
        double espect_escape_private[][nbins*nt_out], particles *ptl_time)
{
    int offset;
    /* Check if the particle already exit the box. */
    if (*iescape_pre == 0 || bc_flag == 0) {
        // advance(&ptl[i], dt);
        tracking_wirz(&ptl[iptl], dt);
        particle_bc(&ptl[iptl].x, &ptl[iptl].y, &ptl[iptl].z, iescape_after);
        nsteps_ptl_tracking[iptl] += 1;
        /* 
         * Trajectory diagnostics. nsteps_output is the time step
         * interval for trajectory info output.
         */
        if (traj_diagnose == 1) {
            if (it % nsteps_output == 0) {
                offset = ntraj_offset_local + (*ntraj_diagnostics_points);
                /* Save the diagnostic particles in aligned memory. */
                memcpy(&ptl_time[offset], &ptl[iptl], sizeof(particles));
                (*ntraj_diagnostics_points)++;
            }
        }
    }
    if (*iescape_pre == 0 && *iescape_after == 1 && bc_flag == 1) {
        /* Only do this when particles escape and open boundary. */
        ptl_energy_fixed(iptl, *ntp, tid, ptl, nbins, nt_out, pmass,
                espect_escape_private);
    }
    int tmp;
    tmp = *iescape_after;
    *iescape_pre = tmp; /* Update particle status. */
    if (it % einterval_t == 0) {
        if (*iescape_after == 0) {
            /* Particles haven't escaped from simulation box. */
            ptl_energy_fixed(iptl, *ntp, tid, ptl, nbins, nt_out, pmass,
                    espect_private);
        }
        (*ntp)++;
    }
}

/******************************************************************************
 * Gather particle energy spectra from different threads. This will calculate
 * particle spectrum inside the box and escaping the simulation box.
 *
 * Input:
 *  it: current time point.
 *  num_threads: total number of threads.
 *  espect_private, espect_escape_private: spectra in each thread.
 * Output:
 *  espectrum, espect_escape: spectra summed over all threads.
 ******************************************************************************/
void GatherParticleSpectra(int num_threads, int nbins, int nt_out,
        double espectrum[][nbins], double espect_private[][nbins*nt_out],
        double espect_escape[][nbins], double espect_escape_private[][nbins*nt_out])
{
    int it, j, k;
    #pragma omp for private(j,k)
    for (it=0; it<nt_out; it++) {
        for (j=0; j<nbins; j++) {
            for (k=0; k<num_threads; k++) {
                espectrum[it][j] += espect_private[k][it*nbins+j];
                if (bc_flag == 1) {
                    /* Only for open boundary conditions. */
                    espect_escape[it][j] += 
                        espect_escape_private[k][it*nbins+j];
                }
            }
        }
    }
}

/******************************************************************************
 * Intilize the arrays for energy spectrum diagnostics.
 ******************************************************************************/
void init_spectrum(int nbins, int nt_out, double espectrum [][nbins])
{
    int i, j;
    for (i = 0; i < nt_out; i++) {
        for (j = 0; j < nbins; j++) {
            espectrum[i][j] = 0.0;
        }
    }
}

/******************************************************************************
 * Particle tracking procedure that uses adaptive step size.
 *
 * Input:
 *  nptl: particle number for this process.
 *  ptl: structure array for particles.
 * Output:
 *  espectrum, espect_tot, espect_escape, espect_private, espect_escape_private
 *  will be updateds.
 ******************************************************************************/
/* void particle_tracking_adaptive(int nptl, */ 
/*         int traj_diagnose, struct particles *ptl, */ 
/*         double espectrum[][nbins], double espect_tot[][nbins], */ 
/*         double espect_escape[][nbins], double espect_private[][nbins*nt_out], */ 
/*         double espect_escape_private[][nbins*nt_out]) */
/* { */
/*     double tps[nt_out]; */
/*     double tott = simul_domain.tmax-simul_domain.tmin; */
/*     time_points(tott, 0, nt_out, tps); */

/*     int ptl_id; */
/*     /1* Number of steps each particle is tracked. *1/ */
/*     for (ptl_id=0; ptl_id<nptl; ptl_id++) { */
/*         nsteps_ptl_tracking[ptl_id] = 1; */
/*     } */

/*     int num_threads = omp_get_max_threads(); */

/*     #pragma omp parallel */
/*     { */
/*         int i, j; */
/*         int tid = omp_get_thread_num(); // current thread ID. */
/*         if (tid == 0) { */
/*             printf("Number of threads: %d\n", omp_get_num_threads()); */
/*         } */
/*         for (j=0; j<nbins*nt_out; j++) { */
/*             espect_private[tid][j] = 0.0; */
/*             espect_escape_private[tid][j] = 0.0; */
/*         } */
/* //        printf("tid %d\n", tid); */
/*         #pragma omp for private(j) */
/*         for (i=0; i<nptl; i++) { */
/*             int iescape_pre = 0; */
/*             int iescape_after = 0; */
/*             int tp_pre, tp_after = 1; /1* Tracking dense output points. *1/ */
/*             double y[nvar], dydx[nvar], ydense[nvar*nt_out]; */
/*             int ivar; */
/*             for (ivar=0; ivar<nvar; ivar++) { */
/*                 dydx[ivar] = 0.0; */
/*             } */
/*             for (j=0; j<nvar*nt_out; j++) { */
/*                 ydense[j] = 0.0; */
/*             } */
/*             double t, dtnext; */
/*             y[0] = ptl[i].x; y[1] = ptl[i].y; y[2] = ptl[i].z; */
/*             y[3] = ptl[i].vx; y[4] = ptl[i].vy; y[5] = ptl[i].vz; */
/*             t = ptl[i].t; */
/*             derivs(t, y, dydx); */
/*             dtnext = 1.0E-5; */

/*             int ntraj_offset_local = 0, offset; */
/*             if (traj_diagnose == 1) { */
/*                 /1* Trajectory diagnostics. *1/ */
/*                 ntraj_offset_local = 0; */
/*                 if (i != 0) { */
/*                     ntraj_offset_local = ntraj_accum[i-1]; */
/*                 } */
/*             } */

/*             int ntraj_diagnostics_points = 1; */

/*             while (t < simul_domain.tmax) { */
/*                 tp_pre = tp_after; */
/*                 if (iescape_pre == 0 || bc_flag == 0) { */
/*                     /1* Check if the particle already exit the box. */
/*                      * Using Bulirsch-Stoer ODE Integrator. *1/ */
/*                     StepperBS(y, dydx, &t, &dtnext, tps, ydense, &tp_after, true); */
/*                     ptl[i].x = y[0]; ptl[i].y = y[1]; ptl[i].z = y[2]; */
/*                     ptl[i].vx = y[3]; ptl[i].vy = y[4]; ptl[i].vz = y[5]; */
/*                     ptl[i].t = t; */
/*                     particle_bc(&ptl[i].x, &ptl[i].y, &ptl[i].z, &iescape_after); */
/*                     nsteps_ptl_tracking[i]++; */
/*                     /1* */ 
/*                      * Trajectory diagnostics. nsteps_output is the time step */
/*                      * interval for trajectory info output. */
/*                      *1/ */
/*                     if (traj_diagnose == 1) { */
/*                         if (j % nsteps_output == 0) { */
/*                             offset = ntraj_offset_local + ntraj_diagnostics_points; */
/*                             memcpy(&ptl_time[offset], &ptl[i], sizeof(particles)); */
/*                             ntraj_diagnostics_points++; */
/*                         } */
/*                     } */
/*                 } */
/*                 if (iescape_pre == 0 && iescape_after == 1 && bc_flag == 1) { */
/*                     /1* Only do this when particles escape and open boundary. *1/ */
/*                     ptl_energy_adaptive(ydense, tp_after, */ 
/*                             tp_pre, tid, espect_escape_private); */
/*                 } */
/*                 iescape_pre = iescape_after; /1* Update particle status. *1/ */
/*                 if (tp_after > tp_pre && iescape_after == 0) { */
/*                     /1* It passes dense output points, so dense output is needed. *1/ */
/*                     ptl_energy_adaptive(ydense, tp_after, */ 
/*                             tp_pre, tid, espect_private); */
/*                 } */
/*             } */
/*         } */

/*         int k; */
/*         #pragma omp for private(j,k) */
/*         for (i=0; i<nt_out; i++) { */
/*             for (j=0; j<nbins; j++) { */
/*                 for (k=0; k<num_threads; k++) { */
/*                     espectrum[i][j] += espect_private[k][i*nbins+j]; */
/*                     if (bc_flag == 1) { */
/*                         /1* Only for open boundary conditions. *1/ */
/*                         espect_escape[i][j] += espect_escape_private[k][i*nbins+j]; */
/*                     } */
/*                 } */
/*             } */
/*         } */
/*     } */
/* } */


/******************************************************************************
 * Particle tracking procedure that can both fixed step and adaptive step size.
 * The time evolution of particle energy spectra are tracked.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nptl: particle number for this process.
 *  dt: time step for particle tracking.
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *  bc_flag: boundary condition flag. 0 for periodic; 1 for open.
 *  tracking_method: 0 for fixed step. 1 for adaptive step.
 *  traj_diagnose: flag for whether to do trajectory diagnostics. 
 *  simul_domain: the simulation domain information.
 *
 * Output:
 *  espectrum: energy spectrum of particles.
 *  ptl: particles struct array will be updated.
 *
 * Input & output:
 *  nsteps_ptl_tracking: total number of tracking steps for each particle.
 ******************************************************************************/
void particle_tracking_hybrid(int mpi_rank, int nptl, double dt, int nbins,
        int nt_out, int bc_flag, int nsteps_output, double pmass,
        int *ntraj_accum, int tracking_method, int traj_diagnose,
        int *nsteps_ptl_tracking, particles *ptl, particles *ptl_time)
{
    double espectrum[nt_out][nbins], espect_tot[nt_out][nbins];
    /* Energy spectrum for escaping particles. */
    double espect_escape[nt_out][nbins]; 
    init_spectrum(nt_out, nbins, espectrum);
    init_spectrum(nt_out, nbins, espect_tot);
    init_spectrum(nt_out, nbins, espect_escape);

    double t1 = omp_get_wtime();
    int num_threads = omp_get_max_threads();
    double espect_private[num_threads+1][nbins*nt_out];
    double espect_escape_private[num_threads+1][nbins*nt_out];
    if (tracking_method == 0) {
        particle_tracking_fixed(nptl, dt, nbins, nt_out, traj_diagnose,
                nsteps_output, pmass, nsteps_ptl_tracking, ptl, ntraj_accum,
                espectrum, espect_tot, espect_escape, espect_private,
                espect_escape_private, ptl_time);
    } else if (tracking_method == 1) {
        /* particle_tracking_adaptive(nptl, traj_diagnose, ptl, espectrum, */ 
        /*         espect_tot, espect_escape, espect_private, espect_escape_private); */
    }
    double t2 = omp_get_wtime();
    printf("Elapsed time: %lf s\n", t2-t1);
    if (traj_diagnose == 0) {
        collect_espectrum(mpi_rank, nbins, nt_out, espectrum, espect_tot,
                "data/espectrum.dat");
        if (bc_flag == 1) {
            /* Only for open boundary conditions. */
            collect_espectrum(mpi_rank, nbins, nt_out, espectrum, espect_tot,
                    "data/espect_escape.dat");
        }
    }
}

/******************************************************************************
 * Dense output time points for adaptive particle tracking method.
 *
 * Input:
 *  tott: total particle tracking time.
 *  nt: # of time points.
 *  iexpo: whether to exponential time series. 0 for no, 1 for yes.
 *  nt_out: number of diagnostic time frames.
 *
 * Output:
 *  tps: the array for those time points.
 ******************************************************************************/
void time_points(double tott, int iexpo, int nt, int nt_out, double *tps)
{
    const double tmin = 1.0E-7; // Small time
    double dt, dt_log, tmin_log;
    int i;
    dt = tott/nt;
    tmin_log = log(tmin);
    dt_log = (log(tott)-tmin_log)/nt_out;
    if (iexpo == 0) {
        for (i=0; i<nt_out; i++) {
            tps[i] = dt * (i+1);
            //printf("i, tps: %d %f\n", i, tps[i]);
        }
    }
    else if (iexpo == 1) {
        for (i=0; i<nt_out; i++) {
            tps[i] = exp(tmin_log+dt_log*i); 
        }
    } else {
        printf("ERROR: using exponential or linear time series?");
        exit(0);
    }
}

/******************************************************************************
 * ============== Particle Tracking based on Wirz's method ==============
 * Reference:
 * Comparison of Charged Particle Tracking Methods for Non-Uniform Magnetic Fields,
 * Hann-Shin Mao and Richard E. Wirz, June 2011
 *
 * Wirz, R., Discharge plasma processes of ring-cusp ion thrusters, 
 * Ph.D. Thesis, California Institute of Technology, 2005.
 *
 * Author: Xiaocan Li
 * Date: Oct-26-2012
 ******************************************************************************/ 
void tracking_wirz(struct particles *ptl, double dt)
{
    double h, Btot;
    double vx_minus, vx_plus, vxold;
    double vy_minus, vy_plus, vyold;
    double vz_minus, vz_plus, vzold;
    double deltax_m, deltay_m, deltaz_m;
    double deltax_new, deltay_new, deltaz_new;
    double hparax, hparay, hparaz; // prallel unit vector to B
    double hperpx, hperpy, hperpz; // perpendicular unit vector
    double hrx, hry, hrz; // radial unit vector
    double dx_para, dy_para, dz_para, dx_perp, dy_perp, dz_perp, dx_r, dy_r, dz_r;
    double dtheta, gyroR, omegac; // gyro frequency and Larmor radius
    double xpre, ypre, zpre; // predicted midpoint
    double vperp; // perpendicular velocity
    double vdothp; // dot product of v and hp
    double x0, y0, z0, t0;
    double deltax, deltay, deltaz;
    double sindtheta, cosdtheta;
    double gama;
    struct emfields emf;
    
    x0 = ptl->x; y0 = ptl->y; z0 = ptl->z; t0 = ptl->t;
    vxold = ptl->vx; vyold = ptl->vy; vzold = ptl->vz;
    getdelta(&deltax, &deltay, &deltaz, vxold, vyold, vzold);	
    gama = getgamma(deltax, deltay, deltaz);
//
    emf.Bx = 0.0; emf.By = 0.0; emf.Bz = 0.0;
    emf.Ex = 0.0; emf.Ey = 0.0; emf.Ez = 0.0;

    get_emf(x0, y0, z0, t0, &emf);
    Btot = sqrt(emf.Bx*emf.Bx + emf.By*emf.By + emf.Bz*emf.Bz);
    h = dt / 2.0;
//    printf("%19.18e\n", dt);
//
// Calculate the first intermidate velocity v_minus by applying the first electric half impulse
    deltax_m = deltax + charge_mass*h*emf.Ex;
    deltay_m = deltay + charge_mass*h*emf.Ey;
    deltaz_m = deltaz + charge_mass*h*emf.Ez;
//
//    printf("deltax1 deltay1 deltaz1 %20.19e %20.19e %20.19e\n ", deltax, deltay, deltaz);
//    printf("deltax_m deltay_m deltaz_m %20.19e %20.19e %20.19e\n ", deltax_m, deltay_m, deltaz_m);
//
    gama = getgamma(deltax, deltay, deltaz);
    vxold = deltax / gama;
    vyold = deltay / gama;
    vzold = deltaz / gama;

    gama = getgamma(deltax_m, deltay_m, deltaz_m);
    vx_minus = deltax_m / gama;
    vy_minus = deltay_m / gama;
    vz_minus = deltaz_m / gama;

// vminus is used to calculate a predicted midpoint
    hparax = emf.Bx / Btot;
    hparay = emf.By / Btot;
    hparaz = emf.Bz / Btot;
    vdothp = vx_minus*hparax + vy_minus*hparay + vz_minus*hparaz; 
    vperp = modulus(vx_minus-vdothp*hparax, vy_minus-vdothp*hparay, 
                    vz_minus-vdothp*hparaz);
    hperpx = (vx_minus-vdothp*hparax) / vperp;  
    hperpy = (vy_minus-vdothp*hparay) / vperp;  
    hperpz = (vz_minus-vdothp*hparaz) / vperp; 

    cross_product(hparax, hparay, hparaz, hperpx, hperpy, hperpz, &hrx, &hry, &hrz);
    hrx *= charge_sign;
    hry *= charge_sign;
    hrz *= charge_sign;
    
    gyroR = (1.0/charge_mass)*vperp*gama*c0 / Btot / L0; // Larmor radius
    omegac = charge_mass*Btot / gama; // gyro frequency
    dtheta = omegac * dt;
    dx_para = h * vdothp * c0 * hparax / L0;
    dy_para = h * vdothp * c0 * hparay / L0;
    dz_para = h * vdothp * c0 * hparaz / L0;

    sindtheta = sin(dtheta*0.5);
    cosdtheta = cos(dtheta*0.5);
    
    dx_perp = gyroR * sindtheta * hperpx;
    dy_perp = gyroR * sindtheta * hperpy;
    dz_perp = gyroR * sindtheta * hperpz;
    
    dx_r = -gyroR * (1.0-cosdtheta) * hrx;
    dy_r = -gyroR * (1.0-cosdtheta) * hry;
    dz_r = -gyroR * (1.0-cosdtheta) * hrz;

// predicted midpoint
    xpre = x0 + dx_para + dx_perp + dx_r;
    ypre = y0 + dy_para + dy_perp + dy_r;
    zpre = z0 + dz_para + dz_perp + dz_r;
    
// updata the coodinate using the predicted midpoint
    get_emf(xpre, ypre, zpre, t0+h, &emf);
    Btot = sqrt(emf.Bx*emf.Bx+emf.By*emf.By+emf.Bz*emf.Bz);

// vminus is used to calculate a predicted midpoint
    hparax = emf.Bx / Btot;
    hparay = emf.By / Btot;
    hparaz = emf.Bz / Btot;
    
    vdothp = vx_minus*hparax + vy_minus*hparay + vz_minus*hparaz; 
    vperp = modulus(vx_minus-vdothp*hparax, vy_minus-vdothp*hparay, 
                    vz_minus-vdothp*hparaz);
    hperpx = (vx_minus-vdothp*hparax) / vperp;  
    hperpy = (vy_minus-vdothp*hparay) / vperp;  
    hperpz = (vz_minus-vdothp*hparaz) / vperp; 

    cross_product(hparax, hparay, hparaz, hperpx, hperpy, hperpz, &hrx, &hry, &hrz);
    hrx *= charge_sign;
    hry *= charge_sign;
    hrz *= charge_sign;
//    printf("hr %10.9e\n ", hperpx * hrx +  hperpy * hry + hperpz * hrz);
    omegac = charge_mass * Btot / gama; // gyro frequency
    dtheta = omegac * dt;
    
    sindtheta = sin(dtheta);
    cosdtheta = cos(dtheta);

    vx_plus = hparax*vdothp + vperp*(cosdtheta*hperpx-sindtheta*hrx);
    vy_plus = hparay*vdothp + vperp*(cosdtheta*hperpy-sindtheta*hry);
    vz_plus = hparaz*vdothp + vperp*(cosdtheta*hperpz-sindtheta*hrz);
//
// update new velocity
// 
    getdelta(&deltax, &deltay, &deltaz, vx_plus, vy_plus, vz_plus);	
    deltax_new = deltax + charge_mass * h * emf.Ex;
    deltay_new = deltay + charge_mass * h * emf.Ey;
    deltaz_new = deltaz + charge_mass * h * emf.Ez;
    
    gama = getgamma(deltax_new, deltay_new, deltaz_new);
    ptl->vx = deltax_new / gama;
    ptl->vy = deltay_new / gama;
    ptl->vz = deltaz_new / gama;
    
//    printf("%20.19e %20.19e %20.19e\n ", vxold * vxold + vyold * vyold + vzold *vzold,
//            vx_minus * vx_minus + vy_minus * vy_minus + vz_minus * vz_minus,
//            *vx * *vx + *vy * *vy + *vz * *vz);
//
// update new position
//
    vdothp = vx_plus*hparax + vy_plus*hparay + vz_plus*hparaz; 
    vperp = modulus(vx_plus-vdothp*hparax, vy_plus-vdothp*hparay, 
                    vz_plus-vdothp*hparaz);
    gyroR = 1.0 / charge_mass * vperp * gama * c0 / L0 / Btot; // Larmor radius
    dx_para = dt * vdothp * c0 * hparax / L0;
    dy_para = dt * vdothp * c0 * hparay / L0;
    dz_para = dt * vdothp * c0 * hparaz / L0;
    
    dx_perp = gyroR * sindtheta * hperpx;
    dy_perp = gyroR * sindtheta * hperpy;
    dz_perp = gyroR * sindtheta * hperpz;
    
    dx_r = -gyroR * (1.0-cosdtheta) * hrx;
    dy_r = -gyroR * (1.0-cosdtheta) * hry;
    dz_r = -gyroR * (1.0-cosdtheta) * hrz;
    
    ptl->x += dx_para + dx_perp + dx_r;
    ptl->y += dy_para + dy_perp + dy_r;
    ptl->z += dz_para + dz_perp + dz_r;
    ptl->t += dt;
}

/******************************************************************************
 * 4th order Runge-Kutta method to track particles.
 ******************************************************************************/
void advance(struct particles *ptl, double dt)
{
    double h;
    struct emfields emf;
    double vx0, vx1, vx2, vx3;
    double vy0, vy1, vy2, vy3;
    double vz0, vz1, vz2, vz3;
    double x0, x1, x2, x3;
    double y0, y1, y2, y3;
    double z0, z1, z2, z3;
    double t0, t1, t2, t3;
    double deltax, deltay, deltaz;
    double deltax1, deltay1, deltaz1;
    double deltax2, deltay2, deltaz2;
    double deltax3, deltay3, deltaz3;
    double Cx, Cy, Cz; // cross product in Lorentz force
    double gama;
    double RKx1, RKx2, RKx3, RKx4;
    double RKy1, RKy2, RKy3, RKy4;
    double RKz1, RKz2, RKz3, RKz4;
    double RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4;
    double RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4;
    double RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4;
    
    //double B_back; // background B-field
    
    x0 = ptl->x; y0 = ptl->y; z0 = ptl->z;
    vx0 = ptl->vx; vy0 = ptl->vy; vz0 = ptl->vz;
    t0 = ptl->t;
    getdelta(&deltax, &deltay, &deltaz, vx0, vy0, vz0);	

    emf.Bx = 0.0; emf.By = 0.0; emf.Bz = 0.0;
    emf.Ex = 0.0; emf.Ey = 0.0; emf.Ez = 0.0;
//
// ===== now RK1 ====
//
    get_emf(x0, y0, z0, t0, &emf);
    h = dt / 2.0;
//    printf("%19.18e\n", dt);

    cross_product(vx0, vy0, vz0, emf.Bx, emf.By, emf.Bz, &Cx, &Cy, &Cz);
    RKx1 = vx0 * c0 / L0;
    RKy1 = vy0 * c0 / L0;
    RKz1 = vz0 * c0 / L0;
    RKdeltax1 = charge_mass*(emf.Ex + Cx);
    RKdeltay1 = charge_mass*(emf.Ey + Cy);
    RKdeltaz1 = charge_mass*(emf.Ez + Cz);
    
//    printf("RK1 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);
//
// ==== now RK2 ====
//
    x1 = x0 + h * RKx1;
    y1 = y0 + h * RKy1;
    z1 = z0 + h * RKz1;
    t1 = t0 + h;
    
    deltax1 = deltax + h * RKdeltax1;
    deltay1 = deltay + h * RKdeltay1;
    deltaz1 = deltaz + h * RKdeltaz1;
    
    gama = getgamma(deltax1, deltay1, deltaz1);
    vx1 = deltax1 / gama;
    vy1 = deltay1 / gama;
    vz1 = deltaz1 / gama;
    

    get_emf(x1, y1, z1, t1, &emf);
    
    //printf("%20.14e %20.14e %20.14e %20.14e\n", x1, y1, z1, t1);

    cross_product(vx1, vy1, vz1, emf.Bx, emf.By, emf.Bz, &Cx, &Cy, &Cz);
    RKx2 = vx1 * c0 / L0;
    RKy2 = vy1 * c0 / L0;
    RKz2 = vz1 * c0 / L0;
    RKdeltax2 = charge_mass*(emf.Ex + Cx);
    RKdeltay2 = charge_mass*(emf.Ey + Cy);
    RKdeltaz2 = charge_mass*(emf.Ez + Cz);
    
//    printf("RK2 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);
// ==== now RK3 ====
    
    x2 = x0 + h*RKx2;
    y2 = y0 + h*RKy2;
    z2 = z0 + h*RKz2;
    t2 = t0 + h;
    
    deltax2 = deltax + h*RKdeltax2;
    deltay2 = deltay + h*RKdeltay2;
    deltaz2 = deltaz + h*RKdeltaz2;
    
    gama = getgamma(deltax2, deltay2, deltaz2);
    vx2 = deltax2 / gama;
    vy2 = deltay2 / gama;
    vz2 = deltaz2 / gama;
    
    get_emf(x2, y2, z2, t2, &emf);
    
    cross_product(vx2, vy2, vz2, emf.Bx, emf.By, emf.Bz, &Cx, &Cy, &Cz);
    RKx3 = vx2 * c0 / L0;
    RKy3 = vy2 * c0 / L0;
    RKz3 = vz2 * c0 / L0;
    RKdeltax3 = charge_mass*(emf.Ex + Cx);
    RKdeltay3 = charge_mass*(emf.Ey + Cy);
    RKdeltaz3 = charge_mass*(emf.Ez + Cz);

//    printf("RK3 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);
// ==== now RK4 ====
//
    x3 = x0 + 2.0*h*RKx3;
    y3 = y0 + 2.0*h*RKy3;
    z3 = z0 + 2.0*h*RKz3;
    t3 = t0 + 2.0*h;
    
    deltax3 = deltax + 2.0*h*RKdeltax3;
    deltay3 = deltay + 2.0*h*RKdeltay3;
    deltaz3 = deltaz + 2.0*h*RKdeltaz3;
    
    gama = getgamma(deltax3, deltay3, deltaz3);
    vx3 = deltax3 / gama;
    vy3 = deltay3 / gama;
    vz3 = deltaz3 / gama;
    
    get_emf(x3, y3, z3, t3, &emf);
    
    cross_product(vx3, vy3, vz3, emf.Bx, emf.By, emf.Bz, &Cx, &Cy, &Cz);
    RKx4 = vx3 * c0 / L0;
    RKy4 = vy3 * c0 / L0;
    RKz4 = vz3 * c0 / L0;
    RKdeltax4 = charge_mass*(emf.Ex +Cx);
    RKdeltay4 = charge_mass*(emf.Ey +Cy);
    RKdeltaz4 = charge_mass*(emf.Ez +Cz);
    
//    printf("RK4 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);

////==== finishing RK method. old (x,y,z,vx,vy,vz,t) -> new (x,y,z,vx,vy,vz,t)
//        printf("RKx1, RKx2, RKx3, RKx4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKx1, RKx2, RKx3, RKx4);
//        printf("RKy1, RKy2, RKy3, RKy4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKy1, RKy2, RKy3, RKy4);
//        printf("RKz1, RKz2, RKz3, RKz4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKz1, RKz2, RKz3, RKz4);
//        printf("RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4);
//        printf("RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4);
//        printf("RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4: %10.9e %10.9e %10.9e %10.9e\n",
//        RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4);
    ptl->x += h/3.0 * ( RKx1 + 2.0*RKx2 + 2.0*RKx3 + RKx4);
    ptl->y += h/3.0 * ( RKy1 + 2.0*RKy2 + 2.0*RKy3 + RKy4);
    ptl->z += h/3.0 * ( RKz1 + 2.0*RKz2 + 2.0*RKz3 + RKz4);
    
    deltax += h/3.0 * (RKdeltax1 + 2.0*RKdeltax2 + 2.0*RKdeltax3 + RKdeltax4);
    deltay += h/3.0 * (RKdeltay1 + 2.0*RKdeltay2 + 2.0*RKdeltay3 + RKdeltay4);
    deltaz += h/3.0 * (RKdeltaz1 + 2.0*RKdeltaz2 + 2.0*RKdeltaz3 + RKdeltaz4);
    
    gama = getgamma(deltax, deltay, deltaz);
    
    ptl->vx = deltax / gama;
    ptl->vy = deltay / gama;
    ptl->vz = deltaz / gama;
    ptl->t += dt;
    
}

/******************************************************************************
 * Calculate gama * beta. gama is Lorentz factor. Beta is the ratio of velocity
 * over speed of light.
 ******************************************************************************/
void getdelta(double *delta_x, double *delta_y, double *delta_z,
		double beta_x, double beta_y, double beta_z)
{
    double beta_mag, gama;
    beta_mag = sqrt(beta_x * beta_x + beta_y * beta_y + beta_z * beta_z);
    gama = 1.0/sqrt(1.0 - beta_mag * beta_mag);
//    printf("gama %20.19e\n", gama);
    *delta_x = gama * beta_x;
    *delta_y = gama * beta_y;
    *delta_z = gama * beta_z;
}

/******************************************************************************
 * Calculate gama from gama * beta.
 ******************************************************************************/
double getgamma(double dx, double dy, double dz)
{
    double delta, gamma;
    
    delta = sqrt(dx*dx + dy*dy + dz*dz);
    gamma = sqrt(1.0 + delta*delta);
//    printf("%20.19e\n ", gamma);
    return gamma;
}

/******************************************************************************
 * Some functions used in code, including cross product of two vectors and the
 * modulus of a vector.
 ******************************************************************************/
void cross_product(double beta_x, double beta_y, double beta_z,
		double Bx, double By, double Bz,
		double* Cx, double* Cy, double* Cz)
{
    *Cx = beta_y*Bz - beta_z*By;
    *Cy = beta_z*Bx - beta_x*Bz;
    *Cz = beta_x*By - beta_y*Bx;
}

double modulus(double x, double y, double z)
{
    return sqrt(x*x+y*y+z*z);
};

/******************************************************************************
 * This procedure is to declare 2D arrays.
 * From: 
 * http://pleasemakeanote.blogspot.com/2008/06/2d-arrays-in-c-using-malloc.html
 ******************************************************************************/
double** Make2DDoubleArray(int arraySizeX, int arraySizeY) 
{
    double** theArray;
    size_t i;
    //printf("arraySizeY, arraySizeY: %d %d\n", arraySizeX, arraySizeY);
    theArray = (double**) malloc(arraySizeX*sizeof(double*));
    for (i = 0; i < arraySizeX; i++) {
        theArray[i] = (double*) malloc(arraySizeY*sizeof(double));
    }
    return theArray;
}

/******************************************************************************
 * Boundary conditions for particles. The flag bc_flag is used to choose the
 * boundary condition.
 *
 * bc_flag: 0 for periodic boundary condition. 1 for open boundary.
 *
 * Input:
 *  x, y, z: spatial positions of particles.
 *
 * Output:
 *  iescape: flag to see if particles escape the simulation box.
 *           0 for no. 1 for yes.
 ******************************************************************************/
void particle_bc(double *x, double *y, double *z, int *iescape)
{
    double xmax, ymax, zmax;
    double xmin, ymin, zmin;
    double xdim, ydim, zdim;
    xmax = simul_domain.xmax;
    ymax = simul_domain.ymax;
    zmax = simul_domain.zmax;
    xmin = simul_domain.xmin;
    ymin = simul_domain.ymin;
    zmin = simul_domain.zmin;
    xdim = xmax - xmin;
    ydim = ymax - ymin;
    zdim = zmax - zmin;
    if (bc_flag == 0) {
        if (*x > xmax) *x -= xdim;
        if (*y > ymax) *y -= ydim;
        if (*z > zmax) *z -= zdim;
        if (*x < xmin) *x += xdim;
        if (*y < ymin) *y += ydim;
        if (*z < zmin) *z += zdim;
    }
    else if (bc_flag == 1) {
        if ((*x>xmax) || (*y>ymax) || (*z>zmax) ||
            (*x<xmin) || (*y<ymin) || (*z<zmin)) {
            *iescape = 1;
        }
    }
}
