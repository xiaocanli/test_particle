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
#include <mpi.h>
#include "constants.h"
#include "domain.h"
#include "particle_info.h"
#include "diagnostics.h"
#include "wlcs.h"
#include "velocity_field.h"
#include "magnetic_field.h"
#include "emfields.h"
#include "tracking.h"

/******************************************************************************
 * Main program for test particle simulation.
 ******************************************************************************/
int main(int argc, char **argv)
{
    int ierr, mpi_size, mpi_rank;
    char config_file_name[LEN_MAX];
    int ipvd, err, system_type, tracking_method;
    int nptl_tot, bc_flag, charge_sign, nptl_traj_tot;
    double ptemp, pmass, pcharge, charge_mass, vthe;
    int nptl, nbins, nt_out;
    int is_traj_diagnostic, is_single_vel;

    /* ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &ipvd); */
    /* ierr = MPI_Init_thread(0, 0, MPI_THREAD_MULTIPLE, &ipvd); */
    ierr = MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &ipvd);
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);

    snprintf(config_file_name, sizeof(config_file_name), "%s", "init.dat");
    struct domain simul_domain;
    err = read_domain(mpi_rank, config_file_name, &simul_domain,
            &system_type, &tracking_method);
    if (err < 0) {
        ierr = MPI_Finalize();
        return 0;
    }

    is_traj_diagnostic = 0;
    read_particle_info(mpi_rank, config_file_name, &nptl_tot, &ptemp,
            &pmass, &pcharge, &charge_sign, &vthe, &charge_mass, &bc_flag,
            &is_traj_diagnostic, &nptl_traj_tot, &is_single_vel);

    int *nptl_accumulate = (int *)malloc(sizeof(int)*mpi_size);
    for (int i = 0; i < mpi_size; i++) {
        nptl_accumulate[i] = 0;
    }
    assign_particles(mpi_rank, mpi_size, &nptl_tot, &nptl_traj_tot,
            system_type, &nptl, nptl_accumulate);

    struct particles *ptl, *ptl_init;
    ptl = (struct particles*)malloc(sizeof(struct particles)*nptl);
    ptl_init = (struct particles*)malloc(sizeof(struct particles)*nptl);
    initialize_partilces(mpi_rank, mpi_size, config_file_name, simul_domain,
            nptl, vthe, system_type, nptl_accumulate, ptl, is_single_vel);
    memcpy(ptl_init, ptl, sizeof(particles));

    get_spectrum_info(mpi_rank, config_file_name, &nbins, &nt_out);

    /* Number of steps each particle is tracked. */
    int *nsteps_ptl_tracking = (int *)malloc(sizeof(int)*nptl);

    double dt = 2.0E-5;
    calc_energy_spectrum(mpi_rank, nptl, ptl, nbins, nt_out, pmass, bc_flag);
    save_particles_fields(mpi_rank, nptl, ptl, nptl_tot, nptl_accumulate,
            "data/particles_fields_init.h5", system_type);
    set_ptl_params_tracking(bc_flag, charge_mass, charge_sign);

    int nsteps_output = 10;
    int *ntraj_accum;
    int traj_diagnose = 0;
    particles *ptl_time;
    particle_tracking_hybrid(mpi_rank, nptl, dt, nbins, nt_out, bc_flag,
            nsteps_output, pmass, ntraj_accum, tracking_method, traj_diagnose,
            nsteps_ptl_tracking, ptl, ptl_time, ptl_init);
    if (is_traj_diagnostic) {
        trajectory_diagnostics(mpi_rank, mpi_size, nptl, dt, pmass,
                nptl_traj_tot, system_type, ptl, nptl_accumulate,
                nsteps_ptl_tracking, nbins, nt_out, bc_flag, tracking_method,
                ptl_init);
    }

    /* release_memory(ptl); */

    if (system_type == 1) {
        /* wire-loop current system */
        free_config();
    }
    if (system_type == 2) {
        free_vfield();
    }
    if (system_type == 3) {
        free_vfield();
        free_bfield();
    }
    free(nsteps_ptl_tracking);
    free(ptl);
    free(ptl_init);
    free(nptl_accumulate);
    ierr = MPI_Finalize();
    return 0;
}
