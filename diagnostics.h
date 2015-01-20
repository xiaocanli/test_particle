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
#define emin 1.0E-5 // Minimum energy.
#define emax 1.0E2  // Maximum energy.
#define rest_ene_proton 938.272046 // Rest energy of proton.
#define logemin log(emin)
#define logemax log(emax)
#define logde ((logemax-logemin)/(nbins-2.0))

void energy_spectrum(int nptl, int it, struct particles *ptl);
void particle_info(int iptl, int ntraj_shift, 
        int ntraj, struct particles *ptl);
void ptl_energy_adaptive(double *ydense, int it1, int it2, int tid, 
        double espectrum[][nbins*nt_out]);
void ptl_energy_fixed(int ptl_id, int it, int tid, 
        struct particles *ptl, double espectrum[][nbins*nt_out]);
void espectrum_collect(double espectrum[][nbins], 
        double espect_tot[][nbins], char *fname);
void save_particles_fields(int nptl, struct particles *ptl, 
        int ntot, int *naccumulate, char *fname);
void read_particles(int nptl, struct particles *ptl, char *fname);
void read_rerun_flag(void);
void sort_particles_energy(int nptl, struct particles *ptl);
void particle_broadcast(int ntot, int *nptl);
void init_ptl_traj(int nptl, int ntest_ptl, int *ntraj, 
        struct particles *ptl, struct particles *ptl_traj);
void trajectory_diagnostics(int nptl, double dt, struct particles *ptl);
