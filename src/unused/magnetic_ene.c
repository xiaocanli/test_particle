/******************************************************************************
 * Calculate the magnetic energy density using force-free field.
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "wlcs.h"
#include "emfields.h"

int main(int argc, char *argv[])
{
    int i, j, k;
    double lmax, dxyz, x, y, z;
    double bene = 0.0;
    struct emfields *emf = (emfields *)malloc(sizeof(struct emfields));
    int nmax = 256;
    lmax = 20.0;
    dxyz = lmax / nmax;
    int nbins = 100;
    int *b_dist = (int *)malloc(sizeof(int)*nbins);
    double *barray = (double *)malloc(sizeof(double)*nbins);
    double bmin, bmax, db;
    bmin = 0.0;
    bmax = 1.6;
    db = (bmax-bmin) / nbins;
    for (i=0; i<nbins; i++) {
        b_dist[i] = 0;
        barray[i] = bmin + db*i;
    }

    double bamp;
    int ibin;
    for (i=0; i<nmax; i++) {
        x = dxyz * i;
        for (j=0; j<nmax; j++) {
            y = dxyz * j;
            for (k=0; k<nmax; k++) {
                z = dxyz * k;
                get_emf(x, y, z, 0.0, emf);
                bamp = sqrt(emf->Bx*emf->Bx +
                        emf->By*emf->By + emf->Bz*emf->Bz);
                ibin = floor((bamp-bmin)/db);
                b_dist[ibin]++;
                //if (bamp > bmax) {
                //    bmax = bamp;
                //}
                bene += bamp * bamp;
            }
        }
    }
    //printf("Maximum of B field: %lf\n", bmax);
    bene = bene/nmax/nmax/nmax;
    printf("Averaged magnetic field energy density: %lf\n", bene);
    FILE *fp;
    fp = fopen("data/bene.dat", "w");
    fprintf(fp, "%lf\n", bene);
    fclose(fp);

    fp = fopen("data/bdist.dat", "w");
    for (i=0; i<nbins; i++) {
        fprintf(fp, "%lf %d\n", barray[i], b_dist[i]);
    }
    fclose(fp);
    free(emf);
    free(b_dist);
    free(barray);
    return 0;
}
