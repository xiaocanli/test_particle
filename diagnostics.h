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

/* Maximum number of trajectory points for one ptl. */
#ifndef MAX_TP
#define MAX_TP 1.0E6
#endif

void get_spectrum_info(int mpi_rank, char *config_file_name, int *nbins,
        int *nt_out);
void calc_energy_spectrum(int mpi_rank, int nptl, struct particles *ptl,
        int nbins, int nt_out, double pmass, int bc_flag);
void save_particles_fields(int mpi_rank, int nptl, struct particles *ptl,
        int ntot, int *naccumulate, char *fname, int system_type);
/* void particle_info(int iptl, int ntraj_shift, */ 
/*         int ntraj, struct particles *ptl); */
/* void ptl_energy_adaptive(double *ydense, int it1, int it2, int tid, */ 
/*         double espectrum[][nbins*nt_out]); */
/* void ptl_energy_fixed(int ptl_id, int it, int tid, */ 
/*         struct particles *ptl, double espectrum[][nbins*nt_out]); */
/* void espectrum_collect(double espectrum[][nbins], */ 
/*         double espect_tot[][nbins], char *fname); */
/* void sort_particles_energy(int nptl, struct particles *ptl); */
/* void init_ptl_traj(int nptl, int ntest_ptl, int *ntraj, */ 
/*         struct particles *ptl, struct particles *ptl_traj); */
/* void trajectory_diagnostics(int nptl, double dt, struct particles *ptl); */
