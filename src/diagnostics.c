/******************************************************************************
* This file is part of CHAOTICB.
* Copyright (C) <2012-2014> <Xiaocan Li> <xl0009@uah.edu>
*
* CHAOTICB is free software: you can redistribute it and/or modify
* it under the terms of the GNU General Public License as published by
* the Free Software Foundation, either version 3 of the License, or
* (at your option) any later version.

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
#include "hdf5.h"
#include "constants.h"
#include "diagnostics.h"
#include "particle_info.h"
#include "emfields.h"
#include "quick_sort.h"
#include "tracking.h"

double emin, emax, logemin, logemax, logde;
/* /1* All of the particle trajectories are output in nsteps_output frames *1/ */
/* /1* init_ptl_traj may update it depending on the maximum particle tracking steps *1/ */
/* int nsteps_output = 10; */


/******************************************************************************
 * Calculate particle energy spectrum for the initial condition.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nptl: particle number for this process.
 *  ptl: structure array of particles.
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *  pmass: particle mass in the unit of proton mass.
 *  bc_flag: boundary condition flag.
 ******************************************************************************/
void calc_energy_spectrum(int mpi_rank, int nptl, struct particles *ptl,
        int nbins, int nt_out, double pmass, int bc_flag)
{
    double gama, ene;
    double *ebins, *einterval;
    double *espectrum, *espect_tot;
    int i, ibin;
    FILE *fp;
    espectrum = (double *)malloc(sizeof(double)*nbins);
    for (i = 0; i < nbins; i++) {
        espectrum[i] = 0;
    }
    if (mpi_rank == 0) {
        espect_tot = (double *)malloc(sizeof(double)*nbins);
        ebins = (double *)malloc(sizeof(double)*nbins);
        einterval = (double *)malloc(sizeof(double)*nbins);
        for (i = 0; i < nbins; i++) {
            espect_tot[i] = 0.0;
        }
        ebins[0] = 0.5 * emin;
        einterval[0] = emin;
        for (i = 1; i < nbins-1; i++) {
            ebins[i] = pow(10, logemin+logde*(i-0.5));
            einterval[i] = pow(10, logemin)*(pow(10, logde*i)-
                    pow(10, logde*(i-1)));
        }
        ebins[nbins-1] = pow(10, logemax+0.5*logde);
        einterval[nbins-1] = pow(10, logemin)*(pow(10, logde*nbins) -
                pow(10, logde*(nbins-1)));
    }

    for (i = 0; i < nptl; i++) {
        gama = sqrt(1.0 + ptl[i].ux*ptl[i].ux +
                ptl[i].uy*ptl[i].uy + ptl[i].uz*ptl[i].uz);
        ene = (gama - 1.0);
        ibin = (log10(ene) - logemin)/logde;
        if (ibin < 0) {
            ibin = 0;
        } else if (ibin > nbins-1) {
            ibin = nbins-1;
        }
        espectrum[ibin]++;
    }

    /* Communication with different processes. */
    MPI_Reduce(espectrum, espect_tot, nbins, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (mpi_rank == 0 ) {
        fp = fopen("data/espectrum.dat", "w");
        for (i = 0; i < nbins; i++) {
            fprintf(fp, "%10.4e ", ebins[i]);
        }
        fprintf(fp, "\n");
        for (i = 0; i < nbins; i++) {
            espectrum[i] = espect_tot[i] / einterval[i];
            fprintf(fp, "%10.4e ", espectrum[i]);
        }
        fprintf(fp, "\n");
        fclose(fp);
        /* If it is open boundary condition, tracking energy spectrum */
        /* for escaping particles. */
        if (bc_flag == 1) {
            fp = fopen("data/espect_escape.dat", "w");
            for (i = 0; i < nbins; i++) {
                fprintf(fp, "%10.4e ", ebins[i]);
            }
            fprintf(fp, "\n");
            for (i = 0; i < nbins; i++) {
                fprintf(fp, "%10.4e ", 0.0);
            }
            fprintf(fp, "\n");
            fclose(fp);
        }
        free(espect_tot);
        free(ebins);
        free(einterval);
        printf("Finished accumulating particle energy spectrum.\n");
    }
    free(espectrum);
}

/******************************************************************************
 * Calculate the energy of a single particle, and put it into the spectra array.
 * This function is for adaptive step size simulations.
 *
 * Input:
 *  ydense: dense outputs from particle tracking procedure, containing particle
 *      spatial positions, velocities.
 *  t1, t2: dense output starting and ending time points.
 *  tid: current OpenMP thread ID.
 *  nbins: number of energy bins.
 *  nt_out: number of time output frames.
 *  pmass: particle mass in the unit of proton mass.
 *  nvar: variable for dense output.
 *  pmass: particle mass in the unit of proton mass.
 *
 * Output:
 *  espectrum is updated.
 ******************************************************************************/
void ptl_energy_adaptive(double *ydense, int t1, int t2, int tid, int nbins,
        int nt_out, double pmass, int nvar, double espectrum[][nbins*nt_out])
{
    double beta, gama, ene;
    int i, j, ibin;

    for (i = t2; i < t1; i++) {
        j = nvar * (i-1);
        beta = sqrt(ydense[j+3]*ydense[j+3]+ydense[j+4]*ydense[j+4]+
                ydense[j+5]*ydense[j+5]);
        gama = 1.0 / sqrt(1.0-beta*beta);
        ene = gama - 1.0;
        ibin = floor((log10(ene)-logemin)/logde);
        if (ibin < 0) {
            ibin = 0;
        } else if (ibin > nbins-1) {
            ibin = nbins-1;
        }
        espectrum[tid][i*nbins+ibin]++;
    }
}

/******************************************************************************
 * Calculate the energy of a single particle, and put it into the spectra array.
 * This function is for fixed step size simulations.
 *
 * Input:
 *  ptl_id: particle id
 *  it: the time point when updating the energy spectrum.
 *  tid: current OpenMP thread ID.
 *  ptl: structure array of particles.
 *  nbins: number of energy bins.
 *  nt_out: number of time output frames.
 *  pmass: particle mass in the unit of proton mass.
 *
 * Output:
 *  espectrum is updated.
 ******************************************************************************/
void ptl_energy_fixed(int ptl_id, int it, int tid, struct particles *ptl,
        int nbins, int nt_out, double espectrum[][nbins*nt_out])
{
    double gama, ene;
    int ibin;

    gama = sqrt(1.0 +
            ptl[ptl_id].ux*ptl[ptl_id].ux +
            ptl[ptl_id].uy*ptl[ptl_id].uy +
            ptl[ptl_id].uz*ptl[ptl_id].uz);
    ene = gama - 1.0;
    ibin = floor((log10(ene)-logemin)/logde);
    if (ibin < 0) {
        ibin = 0;
    } else if (ibin > nbins-1) {
        ibin = nbins-1;
    }
    espectrum[tid][it*nbins+ibin]++;
}

/******************************************************************************
 * Collect particle energy spectrum from each MPI process and output the total
 * energy spectrum.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  espectrum: particle energy spectrum from each MPI process.
 *  nbins: number of energy bins.
 *  nt_out: number of time output frames.
 *  fname: the filename for output.
 *
 * Output:
 *  espect_tot: total particle energy spectrum.
 ******************************************************************************/
void collect_espectrum(int mpi_rank, int nbins, int nt_out,
        double espectrum[][nbins], double espect_tot[][nbins], char *fname)
{
    double einterval[nbins];
    FILE *fp;
    int i, j;
    einterval[0] = emin;
    for (i=1; i<nbins-1; i++) {
        einterval[i] = pow(10, logemin)*(pow(10, logde*i)-pow(10, logde*(i-1)));
    }
    einterval[nbins-1] = pow(10, logemin)*(pow(10, logde*nbins) -
            pow(10, logde*(nbins-1)));

    /* Reduce the spectrum to MPI process 0 */
    MPI_Reduce(espectrum, espect_tot, nbins*nt_out, MPI_DOUBLE, MPI_SUM, 0,
            MPI_COMM_WORLD);
    if (mpi_rank == 0 ) {
        fp = fopen(fname, "a");
        for (i=0; i<nt_out; i++) {
            for (j=0; j<nbins; j++) {
                espectrum[i][j] = espect_tot[i][j] / einterval[j];
                fprintf(fp, "%10.4e ", espectrum[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

/******************************************************************************
 * Save particle information, including positions, velocities, electromagnetic
 * fields at where the particle is for all cases including ADAPTIVE-STEP CASES.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  iptl: particle ID.
 *  ntraj_shift: to decide whether it is the first frame.
 *  ntraj: number trajectories for this output.
 *  ptl: structure array of particles.
 ******************************************************************************/
void save_particle_info(int mpi_rank, int iptl, int ntraj_shift, int ntraj,
        struct particles *ptl)
{
    FILE *fp;
    if (mpi_rank == 0) {
        if (ntraj_shift == 0) {
            fp = fopen("data/ptl_info.bin", "w");
        } else {
            fp = fopen("data/ptl_info.bin", "a");
        }
        fseek(fp, sizeof(struct particles)*(ntraj_shift+ntraj), SEEK_SET);
        fwrite(&ptl[iptl], sizeof(struct particles), 1, fp);
        fclose(fp);
    }
}

/******************************************************************************
 * Save all the particle and fields info into hdf5 file, so that the file can
 * be reused.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nptl: number of particle in one MPI process.
 *  ptl: particle structure array.
 *  ntot: total number of particles/trajectory points for trajectory diagnostics.
 *  naccumulate: accumulated particle # across MPI processes.
 *  fname: filename for the output.
 *  system_type: the type of the system.
 ******************************************************************************/
void save_particles_fields(int mpi_rank, int nptl, struct particles *ptl,
        int ntot, int *naccumulate, char *fname, int system_type)
{
    int i;
    int rank = 1;
    char gname[] = "/particles_fields";
    hid_t file_id, group_id;
    hid_t dset_ptl, dset_emf;
    hid_t filespace, memspace;
    hid_t plist_id;
    hsize_t dimsf[rank];
    hsize_t count[rank], offset[rank];
    hid_t memtype_ptl, filetype_ptl;
    hid_t memtype_emf, filetype_emf;

    /* Initialize data */
    struct emfields *emf_ptl; // emf at particles' position.
    emf_ptl = (struct emfields *)malloc(sizeof(emfields)*nptl);
    for (i=0; i<nptl; i++ ) {
        get_emf(ptl[i].x, ptl[i].y, ptl[i].z, ptl[i].t, &emf_ptl[i]);
    }

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
    file_id = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);
    /* Create a new group */
    group_id = H5Gcreate(file_id, gname, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Create compound datatype */
    create_particles_ctype(&memtype_ptl, &filetype_ptl);
    create_fields_ctype(&memtype_emf, &filetype_emf);

    /* Create the dataspace for the dataset */
    dimsf[0] = ntot;
    filespace = H5Screate_simple(rank, dimsf, NULL);

    /* Create the dataset with default properties and
     * close filespace
     */
    dset_ptl = H5Dcreate(group_id, "particles", filetype_ptl, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    dset_emf = H5Dcreate(group_id, "fields", filetype_emf, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    H5Sclose(filespace);

    /* Count and offset in the memory. */
    count[0] = nptl;
    if (mpi_rank == 0) {
        offset[0] = 0;
    } else {
        offset[0] = naccumulate[mpi_rank-1];
    }

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    memspace = H5Screate_simple(rank, count, NULL);
    filespace = H5Dget_space(dset_ptl);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,
            NULL, count, NULL);
    H5Dwrite(dset_ptl, memtype_ptl, memspace, filespace,
            plist_id, ptl);
    H5Sclose(filespace);
    filespace = H5Dget_space(dset_emf);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset,
            NULL, count, NULL);
    H5Dwrite(dset_emf, memtype_emf, memspace, filespace,
            plist_id, emf_ptl);
    H5Sclose(filespace);

    free(emf_ptl);

    H5Sclose(memspace);
    H5Pclose(plist_id);
    H5Dclose(dset_ptl);
    H5Dclose(dset_emf);
    H5Tclose(memtype_ptl);
    H5Tclose(memtype_emf);
    H5Tclose(filetype_ptl);
    H5Tclose(filetype_emf);
    H5Gclose(group_id);
    H5Fclose(file_id);
}

/******************************************************************************
 * Sort the particles according to their energies.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nptl: number of particles.
 *  ptl: structure array of particles.
 *  pmass: particle mass in the unit of proton mass.
 *  nptl_accumulate: global array for accumulated particle number.
 *
 * Output:
 *  ptl: the sorted particles structure array.
 *  nsteps_ptl_tracking is updated to be consistent with the new ptl array..
 ******************************************************************************/
void sort_particles_energy(int mpi_rank, int nptl, double pmass,
        int *nptl_accumulate, struct particles *ptl, int *nsteps_ptl_tracking)
{
    double gama;
    int i;

    double *ptl_ene = (double *)malloc(sizeof(double)*nptl);
    int *index_ptl = (int *)malloc(sizeof(int)*nptl);
    for (i=0; i<nptl; i++) {
        gama = sqrt(1.0 + ptl[i].ux*ptl[i].ux +
                ptl[i].uy*ptl[i].uy + ptl[i].uz*ptl[i].uz);
        ptl_ene[i] = gama - 1.0;
        index_ptl[i] = i;
    }
    quicksort(ptl_ene, index_ptl, nptl);
    struct particles *ptl_tmp =
        (struct particles *)malloc(sizeof(particles)*nptl);

    /* Reload particles from initial file. */
    read_particles(mpi_rank, nptl, nptl_accumulate, ptl_tmp,
            "data/particles_fields_init.h5");

    int *nsteps = (int *)malloc(sizeof(int)*nptl);
    memcpy(nsteps, nsteps_ptl_tracking, sizeof(int)*nptl);
    for (i=0; i<nptl; i++) {
        memcpy(&ptl[i], &ptl_tmp[index_ptl[i]], sizeof(particles));
        nsteps_ptl_tracking[i] = nsteps[index_ptl[i]];
    }
    free(nsteps);
    free(index_ptl);
    free(ptl_ene);
    free(ptl_tmp);
}

/******************************************************************************
 * Choose particles for trajectory diagnostics from the original sorted
 * particle array. The sampled particles are uniformly distributed in the
 * original particle array. The total # of trajectory points of those particles
 * is returned.
 *
 * Input:
 *  nptl: total # of particles in this MPI process.
 *  nptl_traj: # of trajectory diagnostics test particles for current MPI process.
 *  nsteps_ptl_tracking: the number of tracking steps for each particle.
 *  ptl: particle structure array.
 *
 * Output:
 *  ntraj: the total trajectory points for current MPI process.
 *  ntraj_accum: the accumulated trajectory points starting from mpi_rank = 0
 *  ptl_traj: particle structure array for trajectory diagnostics.
 ******************************************************************************/
void init_ptl_traj(int nptl, int nptl_traj, int *nsteps_ptl_tracking,
        particles *ptl, int *nsteps_output, int *ntraj, int *ntraj_accum,
        particles *ptl_traj)
{
    int i, j, interval, ntraj_max;
    interval = nptl / nptl_traj;
    j = nptl - 1; /* Starting from the highest energy */
    ntraj_max = 0;
    for (i = 0; i < nptl_traj; i++) {
        if (nsteps_ptl_tracking[i] > ntraj_max) {
            ntraj_max = nsteps_ptl_tracking[i];
        }
    }
    if (ntraj_max > MAX_TP) {
        *nsteps_output = (*nsteps_output) * ((double)ntraj_max)/MAX_TP;
    }

    *ntraj = 0;
    j = nptl - 1; /* Starting from the highest energy */
    for (i = 0; i < nptl_traj; i++) {
        ntraj_accum[i] = (nsteps_ptl_tracking[j] - 1) / (*nsteps_output) + 1;
        (*ntraj) += ntraj_accum[i];
        memcpy(&ptl_traj[i], &ptl[j], sizeof(particles));
        j -= interval;
    }
    for (i = 1; i < nptl_traj; i++) {
        ntraj_accum[i] += ntraj_accum[i-1];
    }
}

/******************************************************************************
 * Main function for particle trajectory diagnostics.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  dt: time step for fixed time case.
 *  nptl: total # of particles in this MPI process.
 *  pmass: particle mass in the unit of proton mass.
 *  nptl_traj_tot: the total number of particles for trajectory diagnostics.
 *  system_type: the type of the system. (wlcs, forcefree, MHD+test particle)
 *  ptl: particle structure array.
 *  nptl_accumulate: the accumulative particles numbers along mpi_rank.
 *  nsteps_ptl_tracking: the number of tracking steps for each particle.
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *  nsteps_output: time step interval for trajectory diagnostics. It may be
 *                 updated if total number of trajectories points are too large.
 *  bc_flag: boundary condition flag. 0 for periodic; 1 for open.
 *  tracking_method: 0 for fixed step. 1 for adaptive step.
 ******************************************************************************/
void trajectory_diagnostics(int mpi_rank, int mpi_size, int nptl, double dt,
        double pmass, int nptl_traj_tot, int system_type, struct particles *ptl,
        int *nptl_accumulate, int *nsteps_ptl_tracking, int nbins, int nt_out,
        int *nsteps_output, int bc_flag, int tracking_method, particles *ptl_init)
{
    int i, nptl_traj, ntraj;
    int ntraj_offset_local;
    sort_particles_energy(mpi_rank, nptl, pmass, nptl_accumulate, ptl,
            nsteps_ptl_tracking);
    int *nptl_traj_accumulate = (int *)malloc(sizeof(int)*mpi_size);
    for (int i = 0; i < mpi_size; i++) {
        nptl_traj_accumulate[i] = 0;
    }
    particle_broadcast(mpi_rank, mpi_size, nptl_traj_tot, &nptl_traj,
            nptl_traj_accumulate);
    // Particle trajectory information
    struct particles *ptl_traj =
        (struct particles *)malloc(sizeof(particles)*nptl_traj);
    // Number of trajectory points for all particles in current MPI process
    int *ntraj_accum = (int *)malloc(sizeof(int)*nptl_traj);
    init_ptl_traj(nptl, nptl_traj, nsteps_ptl_tracking, ptl, nsteps_output,
            &ntraj, ntraj_accum, ptl_traj);

    // The accumulation of particle trajectories in all MPI processes
    int *ntraj_accum_global = (int *)malloc(sizeof(int)*mpi_size);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&ntraj, 1, MPI_INT, ntraj_accum_global, 1, MPI_INT,
            0, MPI_COMM_WORLD);
    if (mpi_rank == 0) {
        for (i = 1; i < mpi_size; i++) {
            ntraj_accum_global[i] += ntraj_accum_global[i-1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ntraj_accum_global, mpi_size, MPI_INT, 0, MPI_COMM_WORLD);

    // Time evolution of particles.
    particles *ptl_time = (struct particles *)malloc(sizeof(particles)*ntraj);

    // Initial particles info.
    for (i = 0; i < nptl_traj; i++) {
        if (i == 0) {
            ntraj_offset_local = 0;
        } else {
            ntraj_offset_local = ntraj_accum[i-1];
        }
        memcpy(&ptl_time[ntraj_offset_local], &ptl_traj[i],
                sizeof(particles));
    }

    int traj_diagnose = 1;
    particle_tracking_hybrid(mpi_rank, nptl_traj, dt, nbins, nt_out, bc_flag,
            *nsteps_output, pmass, ntraj_accum, tracking_method, traj_diagnose,
            ptl_init, nsteps_ptl_tracking, ptl_traj, ptl_time);

    save_particles_fields(mpi_rank, ntraj, ptl_time,
            ntraj_accum_global[mpi_size-1], ntraj_accum_global,
            "data/particle_diagnostics.h5", system_type);

    free(ptl_time);
    free(ntraj_accum_global);
    free(ntraj_accum);
    free(ptl_traj);
    free(nptl_traj_accumulate);
}

/******************************************************************************
 * Read energy spectrum information.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  config_file_name: the configuration filename.
 *
 * Output:
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *
 * Other spectrum information:
 *  emin: minimum energy in MeV.
 *  emax: maximum energy in MeV.
 *  logemin: log10(emin).
 *  logemax: log10(emax).
 *  logde: the logarithmic scale of energy interval.
 *****************************************************************************/
void get_spectrum_info(int mpi_rank, char *config_file_name, int *nbins,
        int *nt_out)
{
    FILE *fp;
    int msg;
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    fp = fopen(config_file_name, "r");

    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Energy spectrum info") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "nbins: %d\n", nbins);
    if (msg != 1) {
        printf("Failed to read nbins.\n");
        exit(1);
    }

    msg = fscanf(fp, "emin: %lf\n", &emin);
    if (msg != 1) {
        printf("Failed to read emin.\n");
        exit(1);
    }

    msg = fscanf(fp, "emax: %lf\n", &emax);
    if (msg != 1) {
        printf("Failed to read emax.\n");
        exit(1);
    }

    msg = fscanf(fp, "Number of diagnostic time frames: %d\n", nt_out);
    if (msg != 1) {
        printf("Failed to read the number of diagnostic time frames.\n");
        exit(1);
    }

    free(buff);
    fclose(fp);

    logemin = log10(emin);
    logemax = log10(emax);
    logde = (logemax - logemin) / (*nbins - 2);

    if (mpi_rank == 0) {
        printf("============= Spectrum Info ==============\n");
        printf("nbins = %d\n", *nbins);
        printf("Number of diagnostic frames: %d\n", *nt_out);
        printf("emin, emax (in Lorentz factor) = %f %f\n", emin, emax);
        printf("=========================================\n");
    }
}

/******************************************************************************
 * Calculate diffusion coefficients.
 *
 * Input:
 *  ptl_init: the initial particle information.
 *  ptl: current particle information.
 *
 * Output:
 *  drr: spatial diffusion coefficient
 *  dxx: spatial diffusion coefficient along x-direction
 *  dyy: spatial diffusion coefficient along y-direction
 *  dzz: spatial diffusion coefficient along z-direction
 *  dpp: momentum diffusion coefficient
 *  duu: cosine of pitch angle diffusion coefficient
 *  daa: pitch angle diffusion coefficient
 *  dee: energy diffusion coefficient
 *****************************************************************************/
void calc_diff_coeff(particles *ptl_init, particles *ptl, double dt,
        double *drr, double *dxx, double *dyy, double *dzz, double *dpp,
        double *duu, double *daa, double *dee)
{
    double t, vx, vy, vz, vx0, vy0, vz0, gamma, gamma_init, theta, theta_init;
    t = ptl->t * 2 + dt;
    *drr = 0.0;
    *dxx = 0.0;
    *dyy = 0.0;
    *dzz = 0.0;
    *dpp = 0.0;
    *duu = 0.0;
    *daa = 0.0;
    *dee = 0.0;
    *dxx = pow(ptl_init->x - (ptl->x + ptl->xshift), 2) / t;
    *dyy = pow(ptl_init->y - (ptl->y + ptl->yshift), 2) / t;
    *dzz = pow(ptl_init->z - (ptl->z + ptl->zshift), 2) / t;
    *drr = (*dxx) + (*dyy) + (*dzz);
    *dpp = (pow(ptl_init->ux - ptl->ux, 2) + pow(ptl_init->uy - ptl->uy, 2) +
           pow(ptl_init->uz - ptl->uz, 2)) / t;
    gamma = sqrt(1.0 + ptl->ux*ptl->ux + ptl->uy*ptl->uy + ptl->uz*ptl->uz);
    gamma_init = sqrt(1.0 + ptl_init->ux*ptl_init->ux +
            ptl_init->uy*ptl_init->uy + ptl_init->uz*ptl_init->uz);
    *dee = pow(gamma - gamma_init, 2) / t;
    double mu, mu_init, abs_b, abs_u;
    struct emfields emf_ptl, emf_ptl_init;
    get_emf(ptl->x, ptl->y, ptl->z, ptl->t, &emf_ptl);
    get_emf(ptl_init->x, ptl_init->y, ptl_init->z, ptl_init->t, &emf_ptl_init);
    abs_b = sqrt(emf_ptl.Bx*emf_ptl.Bx + emf_ptl.By*emf_ptl.By +
            emf_ptl.Bz*emf_ptl.Bz);
    abs_u = sqrt(ptl->ux*ptl->ux + ptl->uy*ptl->uy + ptl->uz*ptl->uz);
    mu = (emf_ptl.Bx*ptl->ux + emf_ptl.By*ptl->uy +
          emf_ptl.Bz*ptl->uz) / (abs_b * abs_u);
    abs_b = sqrt(emf_ptl_init.Bx*emf_ptl_init.Bx +
                 emf_ptl_init.By*emf_ptl_init.By +
                 emf_ptl_init.Bz*emf_ptl_init.Bz);
    abs_u = sqrt(ptl_init->ux*ptl_init->ux +
                 ptl_init->uy*ptl_init->uy +
                 ptl_init->uz*ptl_init->uz);
    mu_init = (emf_ptl_init.Bx*ptl_init->ux +
               emf_ptl_init.By*ptl_init->uy +
               emf_ptl_init.Bz*ptl_init->uz) / (abs_b * abs_u);
    *duu = pow(mu - mu_init, 2) / t;
    theta = acos(mu) * 180 / M_PI;
    theta_init = acos(mu_init) * 180 / M_PI;
    *daa = pow(theta - theta_init, 2) / t;
}

/******************************************************************************
 * Collect particle diffusion coefficients.
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  nsteps_dcoeffs: the total output time steps for the diffusion coefficients
 *  fname: the filename for output.
 *  drr: spatial diffusion coefficient
 *  dxx: spatial diffusion coefficient along x-direction
 *  dyy: spatial diffusion coefficient along y-direction
 *  dzz: spatial diffusion coefficient along z-direction
 *  dpp: momentum diffusion coefficient
 *  duu: cosine of pitch angle diffusion coefficient
 *  daa: pitch angle diffusion coefficient
 *  dee: energy diffusion coefficient
 ******************************************************************************/
void collect_diff_coeffs(int mpi_rank, int nsteps_dcoeffs, double *drr,
        double *dxx, double *dyy, double *dzz, double *dpp, double *duu,
        double *daa, double *dee, char *fname, int *nptl_remain)
{
    double *drr_global, *dpp_global, *duu_global;
    double *dxx_global, *dyy_global, *dzz_global;
    double *daa_global, *dee_global;
    int *nptl_remain_global;
    FILE *fp;
    int i;
    drr_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    dxx_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    dyy_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    dzz_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    dpp_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    duu_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    daa_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    dee_global = (double *)calloc(nsteps_dcoeffs, sizeof(double));
    nptl_remain_global = (int *)calloc(nsteps_dcoeffs, sizeof(int));
    /* Reduce the spectrum to MPI process 0 */
    MPI_Reduce(drr, drr_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(dxx, dxx_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(dyy, dyy_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(dzz, dzz_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(dpp, dpp_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(duu, duu_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(daa, daa_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(dee, dee_global, nsteps_dcoeffs, MPI_DOUBLE, MPI_SUM,
            0, MPI_COMM_WORLD);
    MPI_Reduce(nptl_remain, nptl_remain_global, nsteps_dcoeffs, MPI_INT,
            MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0 ) {
        fp = fopen(fname, "w");
        fprintf(fp, "drr        dxx        dyy        dzz        dpp");
        fprintf(fp, "        duu        daa        dee        nptl\n");
        for (i=0; i<nsteps_dcoeffs; i++) {
            fprintf(fp, "%10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %10.4e %d\n",
                    drr_global[i], dxx_global[i], dyy_global[i], dzz_global[i],
                    dpp_global[i], duu_global[i], daa_global[i], dee_global[i],
                    nptl_remain_global[i]);
        }
        fclose(fp);
    }
    free(drr_global);
    free(dxx_global);
    free(dyy_global);
    free(dzz_global);
    free(dpp_global);
    free(duu_global);
    free(daa_global);
    free(dee_global);
    free(nptl_remain_global);
}

/******************************************************************************
 * Get time step interval for trajectory diagnostics
 *
 * Input:
 *  mpi_rank: the rank of current MPI process.
 *  config_file_name: the configuration filename.
 *
 * Output:
 *  nbins: number of energy bins.
 *  nt_out: number of diagnostic time frames.
 *****************************************************************************/
void get_tinterval_traj_diagnostics(int mpi_rank, char *config_file_name,
        int *tinterval_traj)
{
    FILE *fp;
    int msg;
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    fp = fopen(config_file_name, "r");

    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Number of particles for trajectory diagnostics:") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "Time step interval for trajectory diagnostics: %d\n", tinterval_traj);
    if (msg != 1) {
        printf("Failed to read time step interval for trajectory diagnostics.\n");
        exit(1);
    }
    free(buff);
    fclose(fp);

    if (mpi_rank == 0) {
        printf("Time step interval for trajectory diagnostics = %d\n", *tinterval_traj);
    }
}
