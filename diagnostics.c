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
#include "hdf5.h"
#include "Global.h"
#include "cbmpi.h"
#include "diagnostics.h"
#include "quick_sort.h"
#include "tracking.h"

/* Maximum number of trajectory points for one ptl. */
#define MAX_TP 1.0E6

int rerun_flag;
int nsteps_output = 10;
int *ntraj_accum;
int *ntraj_accum_global;
struct particles *ptl_time;

void create_fields_ctype(hid_t *memtype, hid_t *filetype);
void create_particles_ctype(hid_t *memtype, hid_t *filetype);
/******************************************************************************
 * Calculate particle energy spectrum for fixed step size.
 *
 * Input:
 *  nptl: particle number for this process.
 *  it: current time point.
 *  ptl: structure array of particles.
 * Output:
 *  espectrum: the particle energy spectrum.
 ******************************************************************************/
void energy_spectrum(int nptl, int it, struct particles *ptl)
{
    //const double de=(emax-emin)/(nbins-2.0); // Energy interval.
    double beta, gama, ene, rest_ene;
    double *ebins, *einterval;
    double *espectrum, *espect_tot;
    int i, ibin;
    FILE *fp;
    rest_ene = rest_ene_proton*pmass;
    espectrum = (double *)malloc(sizeof(double)*nbins);
    for (i=0; i<nbins; i++) {
        espectrum[i] = 0;
    }
    if (mpi_rank == 0) {
        espect_tot = (double *)malloc(sizeof(double)*nbins);
        ebins = (double *)malloc(sizeof(double)*nbins);
        einterval = (double *)malloc(sizeof(double)*nbins);
        for (i=0; i<nbins; i++) {
            espect_tot[i] = 0.0;
        }
        ebins[0] = 0.5*emin;
        einterval[0] = emin;
        for (i=1; i<nbins-1; i++) {
            ebins[i] = exp(logemin+logde*(i-0.5));
            einterval[i] = exp(logemin)*(exp(logde*i)-exp(logde*(i-1)));
        }
        ebins[nbins-1] = exp(logemax+0.5*logde);
        einterval[nbins-1] = exp(logemin)*(exp(logde*nbins)-exp(logde*(nbins-1)));
    }

    for (i=0; i<nptl; i++) {
        beta = sqrt(ptl[i].vx*ptl[i].vx+ptl[i].vy*ptl[i].vy+ptl[i].vz*ptl[i].vz);
        gama = 1.0/sqrt(1.0-beta*beta);
        ene = (gama-1.0)*rest_ene;
        ibin = (log(ene)-logemin)/logde;
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
        if (it == 0) {
            fp = fopen("espectrum.dat", "w");
            for (i=0; i<nbins; i++) {
                fprintf(fp, "%20.14e ", ebins[i]);
            }
            fprintf(fp, "\n");
            for (i=0; i<nbins; i++) {
                espectrum[i] = espect_tot[i] / einterval[i];
                //printf("%20.14e %20.14e\n", einterval[i], espectrum[i]);
                fprintf(fp, "%20.14e ", espectrum[i]);
            }
            fprintf(fp, "\n");
            fclose(fp);
            /* If it is open boundary condition, tracking energy spectrum
             * for escaping particles. */
            if (bc_flag == 1) {
                fp = fopen("espect_escape.dat", "w");
                for (i=0; i<nbins; i++) {
                    fprintf(fp, "%20.14e ", ebins[i]);
                }
                fprintf(fp, "\n");
                for (i=0; i<nbins; i++) {
                    fprintf(fp, "%20.14e ", 0.0);
                }
                fprintf(fp, "\n");
                fclose(fp);
            }
        } 
        //else {
        //    fp = fopen("espectrum.dat", "a");
        //    for (i=0; i<nbins; i++) {
        //        espectrum[i] = espect_tot[i] / einterval[i];
        //        fprintf(fp, "%20.14e ", espectrum[i]);
        //    }
        //    fprintf(fp, "\n");
        //}
        free(espect_tot);
        free(ebins);
        free(einterval);
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
 *  it1, it2: dense output starting and ending time points.
 * Output:
 *  espectrum is updated.
 ******************************************************************************/
void ptl_energy_adaptive(double *ydense, int it1, int it2, int tid, 
        double espectrum[][nbins*nt_out])
{
    double beta, gama, ene, rest_ene;
    int i, j, ibin;

    rest_ene = rest_ene_proton*pmass;
    for (i=it2; i<it1; i++) {
        j = nvar*(i-1);
        beta = sqrt(ydense[j+3]*ydense[j+3]+ydense[j+4]*ydense[j+4]+
                ydense[j+5]*ydense[j+5]);
        gama = 1.0/sqrt(1.0-beta*beta);
        ene = (gama-1.0)*rest_ene;
        ibin = floor((log(ene)-logemin)/logde);
        //printf("Energy bin %d\n", ibin);
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
 * Output:
 *  espectrum is updated.
 ******************************************************************************/
void ptl_energy_fixed(int ptl_id, int it, int tid, struct particles *ptl, 
        double espectrum[][nbins*nt_out])
{
    double beta, gama, ene, rest_ene;
    int ibin;

    rest_ene = rest_ene_proton*pmass;
    beta = sqrt(ptl[ptl_id].vx*ptl[ptl_id].vx + 
            ptl[ptl_id].vy*ptl[ptl_id].vy + 
            ptl[ptl_id].vz*ptl[ptl_id].vz);
    gama = 1.0/sqrt(1.0-beta*beta);
    ene = (gama-1.0)*rest_ene;
    ibin = floor((log(ene)-logemin)/logde);
    //printf("Energy bin %d\n", ibin);
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
 *  espectrum: particle energy spectrum from each MPI process.
 * Output:
 *  espect_tot: total particle energy spectrum.
 ******************************************************************************/
void espectrum_collect(double espectrum[][nbins], 
        double espect_tot[][nbins], char *fname)
{
//    double ebins[nbins];
    double einterval[nbins];
    FILE *fp;
    int i, j;
//    ebins[0] = 0.5*emin;
    einterval[0] = emin;
    for (i=1; i<nbins-1; i++) {
//        ebins[i] = exp(logemin+logde*(i-0.5));
        einterval[i] = exp(logemin)*(exp(logde*i)-exp(logde*(i-1)));
    }
//    ebins[nbins-1] = exp(logemax+0.5*logde);
    einterval[nbins-1] = exp(logemin)*(exp(logde*nbins)-exp(logde*(nbins-1)));

    /* Reduce the spectrum to MPI process 0 */
    MPI_Reduce(espectrum, espect_tot, nbins*nt_out, 
            MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    if (mpi_rank == 0 ) {
        fp = fopen(fname, "a");
//        for (j=0; j<nbins; j++) {
//            fprintf(fp, "%20.14e ", ebins[j]);
//        }
//        fprintf(fp, "\n");
        for (i=0; i<nt_out; i++) {
            for (j=0; j<nbins; j++) {
                espectrum[i][j] = espect_tot[i][j] / einterval[j];
                //printf("%20.14e %20.14e\n", einterval[j], espectrum[i][j]);
                fprintf(fp, "%20.14e ", espectrum[i][j]);
            }
            fprintf(fp, "\n");
        }
        fclose(fp);
    }
}

/******************************************************************************
 * Save particle information, including positions, velocities, electromagnetic
 * fields at where the particle is for all cases including ADAPTIVE-STEP CASES.
 * Input:
 *  iptl: particle ID.
 *  ptl: structure array of particles.
 ******************************************************************************/
void particle_info(int iptl, int ntraj_shift, int ntraj, struct particles *ptl)
{
    FILE *fp;
    if (mpi_rank == 0) {
        if (ntraj_shift == 0) {
            fp = fopen("ptl_info.bin", "w");
        } else {
            fp = fopen("ptl_info.bin", "a");
        }
        fseek(fp, sizeof(struct particles)*(ntraj_shift+ntraj), SEEK_SET);
        fwrite(&ptl[iptl], sizeof(struct particles), 1, fp);
        fclose(fp);
        //fp = fopen("ptl_info.bin", "rb");
        //for (i=0; i<nptl; i++) {
        //    fseek(fp, sizeof(struct particles)*(nptl*ntraj+i), SEEK_SET);
        //    fread(&ptl1, sizeof(struct particles), 1, fp);
        //    printf("%lf %lf %lf\n", ptl1.x, ptl1.y, ptl1.z);
        //}
        //fclose(fp);
    }
}

/******************************************************************************
 * Save all the particle and fields info into hdf5 file, so that the file can
 * be reused.
 *
 * Input:
 *  nptl: number of particle in one MPI process.
 *  ptl: particle structure array
 *  it: current time points. it = 0 indicates initial point.
 *  ntot: total number of particles/trajectory points for 
 *      trajectory diagnostics.
 *  naccumulate: accumulated particle # across MPI processes.
 ******************************************************************************/
void save_particles_fields(int nptl, struct particles *ptl, 
        int ntot, int *naccumulate, char *fname)
{
    int i;
    int rank = 1;
//    char fname[] = "particles_fields.h5";
    char gname[] = "/particles_fields";
    hid_t file_id, group_id;
    hid_t dset_ptl, dset_emf;
    hid_t filespace, memspace;
    hid_t plist_id;
    //herr_t status;
    hsize_t dimsf[rank];
    hsize_t count[rank], offset[rank];
    hid_t memtype_ptl, filetype_ptl;
    hid_t memtype_emf, filetype_emf;

    /* Initialize data */
    struct emfields *emf_ptl; /* emf at particles' position. */
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
 * Create a compound data type containing the fields information for HDF5.
 *
 * Output:
 *  memtype: compound datatype for memory.
 *  filetype: compound datatype for the file.
 ******************************************************************************/
void create_fields_ctype(hid_t *memtype, hid_t *filetype)
{
    //herr_t status;
    /* Create the compound datatype in memory */
    *memtype = H5Tcreate(H5T_COMPOUND, sizeof(emfields));
    H5Tinsert(*memtype, "Bx", 
            HOFFSET(emfields, Bx), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "By", 
            HOFFSET(emfields, By), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Bz", 
            HOFFSET(emfields, Bz), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ex", 
            HOFFSET(emfields, Ex), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ey", 
            HOFFSET(emfields, Ey), H5T_NATIVE_DOUBLE);
    H5Tinsert(*memtype, "Ez", 
            HOFFSET(emfields, Ez), H5T_NATIVE_DOUBLE);
    /*
     * Create the compound datatype for the file.  Because the standard
     * types we are using for the file may have different sizes than
     * the corresponding native types, we must manually calculate the
     * offset of each member.
     */
    *filetype = H5Tcreate(H5T_COMPOUND, 8*6);
    H5Tinsert (*filetype, "Bx", 0, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "By", 8, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Bz", 8*2, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ex", 8*3, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ey", 8*4, H5T_IEEE_F64BE);
    H5Tinsert (*filetype, "Ez", 8*5, H5T_IEEE_F64BE);
}

/******************************************************************************
 * Sort the particles according to their energies.
 *
 * Input:
 *  nptl: number of particles.
 *  ptl: structure array of particles.
 * Output:
 *  ptl: the sorted particles structure array.
 *  nsteps_ptl_tracking is updated to be consistent with the new ptl array..
 ******************************************************************************/
void sort_particles_energy(int nptl, struct particles *ptl)
{
    double beta, gama, rest_ene;
    int i;

    rest_ene = rest_ene_proton*pmass;
    double *ptl_ene = (double *)malloc(sizeof(double)*nptl);
    int *index_ptl = (int *)malloc(sizeof(int)*nptl);
    for (i=0; i<nptl; i++) {
        beta = sqrt(ptl[i].vx*ptl[i].vx+ptl[i].vy*ptl[i].vy+ptl[i].vz*ptl[i].vz);
        gama = 1.0/sqrt(1.0-beta*beta);
        ptl_ene[i] = (gama-1.0)*rest_ene;
        index_ptl[i] = i;
    }
    quicksort(ptl_ene, index_ptl, nptl);
    struct particles *ptl_tmp = 
        (struct particles *)malloc(sizeof(particles)*nptl);

    /* Reload particles from initial file. */
    read_particles(nptl, ptl_tmp, "particles_fields_init.h5");

    int *nsteps = (int *)malloc(sizeof(int)*nptl);
//    memcpy(ptl_tmp, ptl, sizeof(particles)*nptl);
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
 *  ntest_ptl: # of trajectory diagnostics test particles 
 *      for current MPI process.
 *  ptl: particle structure array.  
 * Return:
 *  ntraj: the total trajectory points for current MPI process.
 *  ptl_traj: particle structure array for trajectory diagnostics.
 ******************************************************************************/
void init_ptl_traj(int nptl, int ntest_ptl, int *ntraj, 
        struct particles *ptl, struct particles *ptl_traj)
{
    int i, j, interval, ntraj_max;
    interval = nptl / ntest_ptl;
    j = nptl - 1; /* Starting from the highest energy */
    ntraj_max = 0;
    for (i=0; i<ntest_ptl; i++) {
        if (nsteps_ptl_tracking[j] > ntraj_max) {
            ntraj_max = nsteps_ptl_tracking[j];
        }
    }
    if (ntraj_max > MAX_TP) {
        nsteps_output = nsteps_output*((double)ntraj_max)/MAX_TP;
    }

    *ntraj = 0;
    j = nptl - 1; /* Starting from the highest energy */
    for (i=0; i<ntest_ptl; i++) {
        (*ntraj) += nsteps_ptl_tracking[j] / nsteps_output;
        memcpy(&ptl_traj[i], &ptl[j], sizeof(particles));
        ntraj_accum[i] = nsteps_ptl_tracking[j] / nsteps_output;
        j -= interval;
    }
    for (i=1; i<ntest_ptl; i++) {
        ntraj_accum[i] += ntraj_accum[i-1];
    }
}

/******************************************************************************
 * Main function for particle trajectory diagnostics.
 *
 * Input:
 *  dt: time step for fixed time case.
 *  nptl: total # of particles in this MPI process.
 *  ptl: particle structure array.
 ******************************************************************************/
void trajectory_diagnostics(int nptl, double dt, struct particles *ptl)
{
    int i, ntest_ptl, ntraj;
    int ntraj_offset_local;
    sort_particles_energy(nptl, ptl);
    particle_broadcast(ntest_ptl_tot, &ntest_ptl);
    struct particles *ptl_traj = 
        (struct particles *)malloc(sizeof(particles)*ntest_ptl);
    ntraj_accum = (int *)malloc(sizeof(int)*ntest_ptl);
    init_ptl_traj(nptl, ntest_ptl, &ntraj, ptl, ptl_traj);

    ntraj_accum_global = (int *)malloc(sizeof(int)*mpi_size);
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Gather(&ntraj, 1, MPI_INT, ntraj_accum_global, 
            1, MPI_INT, 0, MPI_COMM_WORLD);
    if (mpi_rank==0) {
        for (i=1; i<mpi_size; i++) {
            ntraj_accum_global[i] += ntraj_accum_global[i-1];
        }
    }
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(ntraj_accum_global, mpi_size, MPI_INT, 0, MPI_COMM_WORLD);

    /* Time evolution of particles. */
    ptl_time = (struct particles *)malloc(sizeof(particles)*ntraj);

    /* Initial particles info. */
    for (i=0; i<ntest_ptl; i++) {
        if (i == 0) {
            ntraj_offset_local = 0;
        } else {
            ntraj_offset_local = ntraj_accum[i-1];
        }
        memcpy(&ptl_time[ntraj_offset_local], &ptl_traj[i], 
                sizeof(particles));
    }

    particle_tracking_hybrid(ntest_ptl, dt, 1, ptl_traj);
    save_particles_fields(ntraj, ptl_time, ntraj_accum_global[mpi_size-1], 
            ntraj_accum_global, "particle_diagnostics.h5");

    free(ntraj_accum_global);
    free(ntraj_accum);
    free(ptl_traj);
    free(ptl_time);
}
