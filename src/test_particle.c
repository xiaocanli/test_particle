/******************************************************************************
 * Test particle simulation
******************************************************************************/
#include <stdlib.h>
#include <string.h>
#include "diagnostics.h"
#include "wlcs.h"
#include "tracking.h"

/******************************************************************************
 * Function to do test-particle simulation
 ******************************************************************************/
void test_particle(int mpi_rank, int mpi_size, char *config_file_name,
        int system_type, int tracking_method, domain simul_domain)
{
    int nptl_tot, bc_flag, nptl_traj_tot;
    double pmass, charge_mass, vthe;
    int nptl, nbins, nt_out;
    int is_traj_diagnostic, is_single_vel;

    // Read particle information
    is_traj_diagnostic = 0;
    read_particle_info(mpi_rank, config_file_name, &nptl_tot, &pmass,
            &vthe, &charge_mass, &bc_flag, &is_traj_diagnostic,
            &nptl_traj_tot, &is_single_vel);

    // Set the particle numbers on each process
    int *nptl_accumulate = (int *)malloc(sizeof(int)*mpi_size);
    for (int i = 0; i < mpi_size; i++) {
        nptl_accumulate[i] = 0;
    }
    assign_particles(mpi_rank, mpi_size, &nptl_tot, &nptl_traj_tot,
            system_type, &nptl, nptl_accumulate);

    // Initialize particles
    struct particles *ptl, *ptl_init;
    ptl = (struct particles*)malloc(sizeof(struct particles)*nptl);
    ptl_init = (struct particles*)malloc(sizeof(struct particles)*nptl);
    initialize_particles(mpi_rank, mpi_size, config_file_name, simul_domain,
            nptl, vthe, charge_mass, system_type, nptl_accumulate, ptl,
            is_single_vel);
    memcpy(ptl_init, ptl, sizeof(particles) * nptl);

    get_spectrum_info(mpi_rank, config_file_name, &nbins, &nt_out);

    // Number of steps each particle is tracked.
    int *nsteps_ptl_tracking = (int *)malloc(sizeof(int)*nptl);

    double dt;
    dt = get_tracking_time_interval(mpi_rank, config_file_name, pmass);
    if (system_type == 2) adjust_dt_normI(&dt);
    calc_energy_spectrum(mpi_rank, nptl, ptl, nbins, nt_out, pmass, bc_flag);
    save_particles_fields(mpi_rank, nptl, ptl, nptl_tot, nptl_accumulate,
            "data/particles_fields_init.h5", system_type);
    set_ptl_params_tracking(bc_flag, charge_mass);

    // tracking particles
    int nsteps_output = 10;
    get_tinterval_traj_diagnostics(mpi_rank, config_file_name, &nsteps_output);
    int *ntraj_accum;
    int traj_diagnose = 0;
    particles *ptl_time;
    particle_tracking_hybrid(mpi_rank, nptl, dt, nbins, nt_out, bc_flag,
            nsteps_output, pmass, ntraj_accum, tracking_method, traj_diagnose,
            ptl_init, nsteps_ptl_tracking, ptl, ptl_time);
    if (is_traj_diagnostic) {
        trajectory_diagnostics(mpi_rank, mpi_size, nptl, dt, pmass,
                nptl_traj_tot, system_type, ptl, nptl_accumulate,
                nsteps_ptl_tracking, nbins, nt_out, &nsteps_output,
                bc_flag, tracking_method, ptl_init);
    }

    free(nsteps_ptl_tracking);
    free(ptl);
    free(ptl_init);
    free(nptl_accumulate);
}
