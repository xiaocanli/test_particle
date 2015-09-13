/******************************************************************************
 * Calculate the magnetic energy density using force-free field.
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* Magnetic fields structure */
typedef struct bfields {
	double Bx, By, Bz;
} bfields;

void getb_ff(double x, double y, double z, double t, struct bfields *bmf);

int main(int argc, char *argv[])
{
    int i, j, k;
    double lmax, dxyz, x, y, z;
    double bene = 0.0;
    struct bfields *bmf = (bfields *)malloc(sizeof(struct bfields));
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
                getb_ff(x, y, z, 0.0, bmf);
                bamp = sqrt(bmf->Bx*bmf->Bx + 
                        bmf->By*bmf->By + bmf->Bz*bmf->Bz);
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
    fp = fopen("bene.dat", "w");
    fprintf(fp, "%lf\n", bene);
    fclose(fp);

    fp = fopen("bdist.dat", "w");
    for (i=0; i<nbins; i++) {
        fprintf(fp, "%lf %d\n", barray[i], b_dist[i]);
    }
    fclose(fp);
    free(bmf);
    free(b_dist);
    free(barray);
    return 0;
}

void getb_ff(double x, double y, double z, double t, struct bfields *bmf)
{
    double B0, lambda0;
    double ffA, ffB, ffC;
    B0 = 1.0;
    lambda0 = 1.0;
    ffA = 1.0;
    ffB = sqrt(2.0/3.0);
    ffC = sqrt(1.0/3.0);

    bmf->Bx = ffA*sin(lambda0*z) + ffC*cos(lambda0*y);
    bmf->By = ffB*sin(lambda0*x) + ffA*cos(lambda0*z);
    bmf->Bz = ffC*sin(lambda0*y) + ffB*cos(lambda0*x);

    bmf->Bx *= B0;
    bmf->By *= B0;
    bmf->Bz *= B0;
}

