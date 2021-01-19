/*----------------------------------------------------------------------
This code gets the Poincaré map of a chaotic (or any) magnetic field.
It uses Cartesian coordinates and a free parameter delta to control
the step size.
    dx/Bx = dy/By = dz/Bz = delta
It requires dl = sqrt(dx^2 + dy^2 + dz^2) < some threshold value.
------------------------------------------------------------------------- */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "constants.h"
#include "domain.h"
#include "emfields.h"
#include <hdf5.h>

void next_location(double *fx, double *fy, double *fz, double dstep);
void field_boundary(double *fx, double *fy, double *fz);
int config_poincare_map(int mpi_rank, char *config_file_name);

struct domain simul_domain;
int map_for_single_point; // whether to get map for a single starting point
int face_cross_section;   // which face chosen as cross section (0: x, 1: y, 2:z)
int nhorizontal;  // number of sampling point along the horizontal direction
int nvertical;    // number of sampling point along the vertical direction
int ncross;       // number of crossing points on the map
double fx0, fy0, fz0;     // starting point if map_for_single_point==1
double min_horizontal, max_horizontal; // min and max for horizontal direction if map_for_single_point==0
double min_vertical, max_vertical;     // min and max for vertical direction if map_for_single_point==0
double step_fieldline, fieldline_length_max;

/******************************************************************************
 * Main function for getting a Poincaré map
 * Arguments:
 *  npoints: number of points on the map
 *  xs, ys, zs: starting points for tracing magnetic field lines
 ******************************************************************************/
int poincare_map(int mpi_rank, int mpi_size, char *config_file_name, domain *sdomain)
{
    if (mpi_rank == 0)
        printf("***** Entering Poincare *****\n");

    memcpy(&simul_domain, sdomain, sizeof(domain));
    config_poincare_map(mpi_rank, config_file_name);

    double *xpos, *ypos, *zpos;
    int npoints, nshift;
    if (map_for_single_point) {
        npoints = 1;
        xpos = (double *)malloc(sizeof(double)*(ncross + 1));
        ypos = (double *)malloc(sizeof(double)*(ncross + 1));
        zpos = (double *)malloc(sizeof(double)*(ncross + 1));
        /* Set to a value that is impossible for a actual Poincare map */
        for (long int i = 0; i < ncross + 1; i++) {
            xpos[i] = simul_domain.xmax * 2;
            ypos[i] = simul_domain.ymax * 2;
            zpos[i] = simul_domain.zmax * 2;
        }
    } else {
        npoints = nhorizontal * nvertical / mpi_size;
        nshift = npoints * mpi_rank;
        int res = (nhorizontal * nvertical) % mpi_size;
        if (res != 0) {
            if (mpi_rank < res) {
                npoints++;
                nshift += mpi_rank;
            } else {
                nshift += res;
            }
        }
        xpos = (double *)malloc(sizeof(double)*(ncross + 1)*npoints);
        ypos = (double *)malloc(sizeof(double)*(ncross + 1)*npoints);
        zpos = (double *)malloc(sizeof(double)*(ncross + 1)*npoints);
        /* Set to a value that is impossible for a actual Poincare map */
        for (long int i = 0; i < (ncross + 1)*npoints; i++) {
            xpos[i] = simul_domain.xmax * 2;
            ypos[i] = simul_domain.ymax * 2;
            zpos[i] = simul_domain.zmax * 2;
        }
    }

    if (map_for_single_point) {
        double point_pre, point_init;
        double cross;
        double pos_init[3], pos[3];
        int icross;
        pos_init[0] = fx0;
        pos_init[1] = fy0;
        pos_init[2] = fz0;
        point_init = pos_init[face_cross_section];
        pos[0] = fx0;
        pos[1] = fy0;
        pos[2] = fz0;
        xpos[0] = fx0;
        ypos[0] = fy0;
        zpos[0] = fz0;
        cross = 1.0;
        icross = 1;
        double fieldline_length = 0.0;
        while(icross <= ncross){
            cross = 1.0;
            while (cross > 0 && fieldline_length < fieldline_length_max) {
                point_pre = pos[face_cross_section];
                next_location(&pos[0], &pos[1], &pos[2], step_fieldline);
                cross = (point_pre - point_init) * (pos[face_cross_section] - point_init);
                fieldline_length += step_fieldline;
            }
            if (fieldline_length < fieldline_length_max) {
                xpos[icross] = pos[0];
                ypos[icross] = pos[1];
                zpos[icross] = pos[2];
                ++icross;
            } else {
                printf("Reached maximum length for field-line tracing before %d crosses\n",
                        icross + 1);
                break;
            }
        }
    } else {
        int ipoint;
#pragma omp parallel for
        for (ipoint = nshift; ipoint < nshift + npoints; ipoint++) {
            double point_pre, point_init;
            double cross, delta_h, delta_v;
            double pos_init[3], pos[3];
            int icross, ih, iv, count, shift;
            int threadid = omp_get_thread_num();
            int num_threads = omp_get_max_threads();
            iv = ipoint / nhorizontal;
            ih = ipoint % nhorizontal;
            if (nhorizontal > 1)
                delta_h = (max_horizontal - min_horizontal) / (nhorizontal - 1);
            else
                delta_h = max_horizontal - min_horizontal;
            if (nvertical > 1)
                delta_v = (max_vertical - min_vertical) / (nvertical - 1);
            else
                delta_v = max_vertical - min_vertical;
            count = 0;
            pos_init[0] = fx0;
            pos_init[1] = fy0;
            pos_init[2] = fz0;
            // Fix pos_init for different faces
            for (int i = 0; i != 3; i++) {
                if (i != face_cross_section) {
                    if (count == 0) {
                        pos_init[i] = delta_h * ih + min_horizontal;
                    } else {
                        pos_init[i] = delta_v * iv + min_vertical;
                    }
                    count++;
                }
            }
            point_init = pos_init[face_cross_section];
            pos[0] = pos_init[0];
            pos[1] = pos_init[1];
            pos[2] = pos_init[2];
            shift = (ipoint - nshift) * (ncross + 1);
            xpos[shift] = pos_init[0];
            ypos[shift] = pos_init[1];
            zpos[shift] = pos_init[2];
            cross = 1.0;
            icross = 1;
            double fieldline_length = 0.0;
            while(icross <= ncross){
                cross = 1.0;
                while (cross > 0 && fieldline_length < fieldline_length_max) {
                    point_pre = pos[face_cross_section];
                    next_location(&pos[0], &pos[1], &pos[2], step_fieldline);
                    cross = (point_pre - point_init) * (pos[face_cross_section] - point_init);
                    fieldline_length += step_fieldline;
                }
                if (fieldline_length < fieldline_length_max) {
                    xpos[shift + icross] = pos[0];
                    ypos[shift + icross] = pos[1];
                    zpos[shift + icross] = pos[2];
                    ++icross;
                } else {
                    printf("Reached maximum length for field-line tracing before %d crosses\n",
                            icross + 1);
                    break;
                }
            }
        }
    }

    // save the data
    const int RANK=2;
    hid_t file_id, group_id, dset_id, plist_id;
    hid_t filespace, memspace;
    hsize_t dset_dims[RANK], dcount[RANK], doffset[RANK];
    MPI_Info fileinfo;
    herr_t status;

    if (map_for_single_point) {
        dset_dims[0] = 1;
        dset_dims[1] = ncross + 1;
        dcount[0] = 1;
        dcount[1] = ncross + 1;
        doffset[0] = 0;
        doffset[1] = 0;
    } else {
        dset_dims[0] = nhorizontal * nvertical;
        dset_dims[1] = ncross + 1;
        dcount[0] = npoints;
        dcount[1] = ncross + 1;
        doffset[0] = nshift;
        doffset[1] = 0;
    }

    MPI_Info_create(&fileinfo);
    MPI_Info_set(fileinfo, "romio_cb_read", "automatic");
    MPI_Info_set(fileinfo, "romio_ds_read", "automatic");

    plist_id = H5Pcreate(H5P_FILE_ACCESS);
    H5Pset_fapl_mpio(plist_id, MPI_COMM_WORLD, fileinfo);
    file_id = H5Fcreate("data/poincare_map.h5", H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
    H5Pclose(plist_id);

    group_id = H5Gcreate(file_id, "/poincare_map", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    filespace = H5Screate_simple(RANK, dset_dims, NULL);
    memspace = H5Screate_simple(RANK, dcount, NULL);
    H5Sselect_hyperslab(filespace, H5S_SELECT_SET, doffset, NULL, dcount, NULL);

    plist_id = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

    dset_id = H5Dcreate(group_id, "xpos", H5T_NATIVE_DOUBLE, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, xpos);
    H5Dclose(dset_id);
    dset_id = H5Dcreate(group_id, "ypos", H5T_NATIVE_DOUBLE, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, ypos);
    H5Dclose(dset_id);
    dset_id = H5Dcreate(group_id, "zpos", H5T_NATIVE_DOUBLE, filespace,
            H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dset_id, H5T_NATIVE_DOUBLE, memspace, filespace, plist_id, zpos);
    H5Dclose(dset_id);

    H5Pclose(plist_id);
    H5Sclose(memspace);
    H5Sclose(filespace);
    H5Gclose(group_id);
    H5Fclose(file_id);

    MPI_Info_free(&fileinfo);

    free(xpos);
    free(ypos);
    free(zpos);

    if (mpi_rank == 0)
        printf("***** Leaving Poincare *****\n");

    return 0;
}

/*======================================================================
Standard RK4 method: http://www.myphysicslab.com/runge_kutta.html
xn+1 = xn + h/6*(a + 2b + 2c + d), where
a = f(tn, xn)
b = f(tn + h/2, xn + h/2*a)
c = f(tn + h/2, xn + h/2*b)
d = f(tn + h, xn + h*c)
========================================================================*/
void next_location(double *fx, double *fy, double *fz, double dstep)
{
    double bx1, by1, bz1;
    double bx2, by2, bz2;
    double bx3, by3, bz3;
    double bx4, by4, bz4;
    double ib, xtmp, ytmp, ztmp;
    struct emfields emf;

    emf.Bx = 0.0; emf.By = 0.0; emf.Bz = 0.0;
    emf.Ex = 0.0; emf.Ey = 0.0; emf.Ez = 0.0;

    field_boundary(fx, fy, fz);
    get_magnetic_field(*fx, *fy, *fz, 0, &emf);
    ib = 1.0 / sqrt(emf.Bx * emf.Bx + emf.By * emf.By + emf.Bz * emf.Bz);
    bx1 = emf.Bx * ib;
    by1 = emf.By * ib;
    bz1 = emf.Bz * ib;

    xtmp = *fx + 0.5 * dstep * bx1;
    ytmp = *fy + 0.5 * dstep * by1;
    ztmp = *fz + 0.5 * dstep * bz1;
    field_boundary(&xtmp, &ytmp, &ztmp);
    get_magnetic_field(xtmp, ytmp, ztmp, 0, &emf);
    ib = 1.0 / sqrt(emf.Bx * emf.Bx + emf.By * emf.By + emf.Bz * emf.Bz);
    bx2 = emf.Bx * ib;
    by2 = emf.By * ib;
    bz2 = emf.Bz * ib;

    xtmp = *fx + 0.5 * dstep * bx2;
    ytmp = *fy + 0.5 * dstep * by2;
    ztmp = *fz + 0.5 * dstep * bz2;
    field_boundary(&xtmp, &ytmp, &ztmp);
    get_magnetic_field(xtmp, ytmp, ztmp, 0, &emf);
    ib = 1.0 / sqrt(emf.Bx * emf.Bx + emf.By * emf.By + emf.Bz * emf.Bz);
    bx3 = emf.Bx * ib;
    by3 = emf.By * ib;
    bz3 = emf.Bz * ib;

    xtmp = *fx + dstep * bx3;
    ytmp = *fy + dstep * by3;
    ztmp = *fz + dstep * bz3;
    field_boundary(&xtmp, &ytmp, &ztmp);
    get_magnetic_field(xtmp, ytmp, ztmp, 0, &emf);
    ib = 1.0 / sqrt(emf.Bx * emf.Bx + emf.By * emf.By + emf.Bz * emf.Bz);
    bx4 = emf.Bx * ib;
    by4 = emf.By * ib;
    bz4 = emf.Bz * ib;

    *fx += (dstep / 6) * (bx1 + 2*bx2 + 2*bx3 + bx4);
    *fy += (dstep / 6) * (by1 + 2*by2 + 2*by3 + by4);
    *fz += (dstep / 6) * (bz1 + 2*bz2 + 2*bz3 + bz4);
}

/******************************************************************************
 * Check if a point crosses the boundaries.
 * We assume periodic boundary here.
 ******************************************************************************/
void field_boundary(double *fx, double *fy, double *fz)
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
    if (*fx > xmax)
        *fx -= xdim;
    if (*fy > ymax)
        *fy -= ydim;
    if (*fz > zmax)
        *fz -= zdim;
    if (*fx < xmin)
        *fx += xdim;
    if (*fy < ymin)
        *fy += ydim;
    if (*fz < zmin)
        *fz += zdim;
}

/******************************************************************************
 * Configuration for Poincare map
 ******************************************************************************/
int config_poincare_map(int mpi_rank, char *config_file_name)
{
    char *buff = (char *)malloc(sizeof(*buff)*LEN_MAX);
    int err, msg;
    FILE *fp;

    fp = fopen(config_file_name, "r");

    while (fgets(buff, LEN_MAX, fp) != NULL) {
        //puts(buff);
        if (strstr(buff, "Poincare map configuration") != NULL) {
            break;
        }
    }
    msg = fscanf(fp, "Which face as the cross section (0: x, 1: y, 2:z):%d\n",
            &face_cross_section);
    if (msg != 1) {
        printf("Failed to get the face for cross section.\n");
        exit(1);
    }
    if (mpi_rank == 0) {
        switch (face_cross_section) {
            case 0:
                printf("The face for cross section is x-face (yz plane).\n");
                break;
            case 1:
                printf("The face for cross section is y-face (xz plane).\n");
                break;
            case 2:
                printf("The face for cross section is z-face (xy plane).\n");
                break;
            default:
                printf("ERROR: wrong face for cross section.");
                return -1;
        }
    }

    msg = fscanf(fp, "Field-line tracking step: %lf\n", &step_fieldline);
    if (msg != 1) {
        printf("Failed to get the step for field-line tracking.\n");
        exit(1);
    }
    if (mpi_rank == 0) {
        printf("The step for field-line tracking: %lf\n", step_fieldline);
    }

    msg = fscanf(fp, "Maximum tracking length: %lf\n", &fieldline_length_max);
    if (msg != 1) {
        printf("Failed to get the maximum length for field-line tracking.\n");
        exit(1);
    }
    if (mpi_rank == 0) {
        printf("The maximum length for field-line tracking: %lf\n", fieldline_length_max);
    }

    msg = fscanf(fp, "Number of crossing points: %d\n", &ncross);
    if (msg != 1) {
        printf("Failed to get the number of crossing points.\n");
        exit(1);
    }
    if (mpi_rank == 0) {
        printf("Number of crossing points: %d\n", ncross);
    }

    map_for_single_point = 1;
    msg = fscanf(fp, "Whether to get the map only for a single point:%d\n",
            &map_for_single_point);
    if (msg != 1) {
        printf("Failed to get flag for whether to get the map for only a single point.\n");
        exit(1);
    }

    msg = fscanf(fp, "x-position of the point: %lf\n", &fx0);
    if (msg != 1) {
        printf("Failed to get the x-position of the point.\n");
        exit(1);
    }
    msg = fscanf(fp, "y-position of the point: %lf\n", &fy0);
    if (msg != 1) {
        printf("Failed to get the y-position of the point.\n");
        exit(1);
    }
    msg = fscanf(fp, "z-position of the point: %lf\n", &fz0);
    if (msg != 1) {
        printf("Failed to get the z-position of the point.\n");
        exit(1);
    }

    if (map_for_single_point) {
        if (mpi_rank == 0) {
            printf("Get the map only for a single point.\n");
            printf("The position of this point is: %lf, %lf, %lf\n", fx0, fy0, fz0);
        }
    } else {
        while (fgets(buff, LEN_MAX, fp) != NULL) {
            if (strstr(buff, "for a 2D plane") != NULL) {
                break;
            }
        }
        msg = fscanf(fp, "min_horizontal: %lf\n", &min_horizontal);
        if (msg != 1) {
            printf("Failed to get min_horizontal.\n");
            exit(1);
        }
        msg = fscanf(fp, "max_horizontal: %lf\n", &max_horizontal);
        if (msg != 1) {
            printf("Failed to get max_horizontal.\n");
            exit(1);
        }
        msg = fscanf(fp, "min_vertical: %lf\n", &min_vertical);
        if (msg != 1) {
            printf("Failed to get min_vertical.\n");
            exit(1);
        }
        msg = fscanf(fp, "max_vertical: %lf\n", &max_vertical);
        if (msg != 1) {
            printf("Failed to get max_vertical.\n");
            exit(1);
        }
        msg = fscanf(fp, "nhorizontal: %d\n", &nhorizontal);
        if (msg != 1) {
            printf("Failed to get nhorizontal.\n");
            exit(1);
        }
        msg = fscanf(fp, "nvertical: %d\n", &nvertical);
        if (msg != 1) {
            printf("Failed to get nvertical.\n");
            exit(1);
        }
        if (mpi_rank == 0) {
            printf("Get the maps for points in a 2D plane.\n");
            printf("Min and max of the horizontal direction are: %lf, %lf\n",
                    min_horizontal, max_horizontal);
            printf("Min and max of the vertical direction are: %lf, %lf\n",
                    min_vertical, max_vertical);
            printf("Number of points along the horizontal direction: %d\n", nhorizontal);
            printf("Number of points along the vertical direction: %d\n", nvertical);
            switch (face_cross_section) {
                case 0:
                    printf("The x-position of the plane: %lf\n", fx0);
                    break;
                case 1:
                    printf("The y-position of the plane: %lf\n", fy0);
                    break;
                case 2:
                    printf("The z-position of the plane: %lf\n", fz0);
                    break;
                default:
                    printf("ERROR: wrong face for cross section.");
                    return -1;
            }
        }
    }
    fclose(fp);

    free(buff);
    return 0;
}
