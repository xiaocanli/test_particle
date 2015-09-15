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
#include "particle_info.h"

/* void time_points(double tott, int iexpo, int nt, double *tps); */
/* double** Make2DDoubleArray(int arraySizeX, int arraySizeY); */
void particle_bc(double *x, double *y, double *z, int *iescape);
/* void init_spectrum(double espectrum [][nbins]); */
/* void particle_tracking_adaptive(int nptl, */ 
/*         int traj_diagnose, struct particles *ptl, */ 
/*         double espectrum[][nbins], double espect_tot[][nbins], */ 
/*         double espect_escape[][nbins], double espect_private[][nbins*nt_out], */ 
/*         double espect_escape_private[][nbins*nt_out]); */
/* void particle_tracking_fixed(int nptl, double dt, */ 
/*         int traj_diagnose, struct particles *ptl, */ 
/*         double espectrum[][nbins], double espect_tot[][nbins], */ 
/*         double espect_escape[][nbins], double espect_private[][nbins*nt_out], */ 
/*         double espect_escape_private[][nbins*nt_out]); */
void tracking_wirz(struct particles *ptl, double dt);

double getgamma(double dx, double dy, double dz);
void getdelta(double *delta_x, double *delta_y, double *delta_z,
        double beta_x, double beta_y, double beta_z);
double modulus(double x, double y, double z);
/* void particle_tracking_hybrid(int nptl, double dt, */ 
/*         int traj_diagnose, struct particles *ptl); */

void set_variables_tracking(grids *sgrid, domain *sdomain, int stype, int mt);
void set_ptl_params_tracking(int pbc, double cmass, int csign);

void cross_product(double beta_x, double beta_y, double beta_z, double Bx,
        double By, double Bz, double* Cx, double* Cy, double* Cz);

void SingleStepParticleTrackingOMP(int iptl, int it, double dt, int tid,
        int einterval_t, int traj_diagnose, int nbins, int nt_out,
        int nsteps_output, double pmass, int *ntp, particles *ptl,
        int *iescape_pre, int *iescape_after, int ntraj_offset_local,
        int *ntraj_diagnostics_points, int *nsteps_ptl_tracking,
        double espect_private[][nbins*nt_out], 
        double espect_escape_private[][nbins*nt_out], particles *ptl_time);

void MultiStepsParticleTrackingOmp(int nptl, double dt, int tid, int sstep,
        int estep, int nbins, int nt_out, int einterval_t, int traj_diagnose,
        int nsteps_output, double pmass, int *ntraj_accum,
        int *ntraj_diagnostics_points_array, struct particles *ptl,
        int *nsteps_ptl_tracking, double espect_private[][nbins*nt_out], 
        double espect_escape_private[][nbins*nt_out], particles *ptl_time);

void particle_tracking_fixed(int nptl, double dt, int nbins, int nt_out,
        int traj_diagnose, int nsteps_output, double pmass,
        int *nsteps_ptl_tracking, particles *ptl, int *ntraj_accum,
        double espectrum[][nbins], double espect_tot[][nbins],
        double espect_escape[][nbins], double espect_private[][nbins*nt_out],
        double espect_escape_private[][nbins*nt_out], particles *ptl_time);

void GatherParticleSpectra(int num_threads, int nbins, int nt_out,
        double espectrum[][nbins], double espect_private[][nbins*nt_out],
        double espect_escape[][nbins], double espect_escape_private[][nbins*nt_out]);

void particle_tracking_hybrid(int mpi_rank, int nptl, double dt, int nbins,
        int nt_out, int bc_flag, int nsteps_output, double pmass,
        int *ntraj_accum, int tracking_method, int traj_diagnose,
        int *nsteps_ptl_tracking, particles *ptl, particles *ptl_time);
