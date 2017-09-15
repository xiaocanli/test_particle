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

/******************************************************************************
 * Bulirsch-Stoer ODE Integrator                                      
 * Reference:                                                          
 * Press, William H., et al. Numerical recipes 3rd edition:            
 * The art of scientific computing. Cambridge University Press, 2007. 
 * Chapter 17                                                         
 * Date:   04-17-2013                                                 
 * Author: Xiaocan Li                                                 
 ******************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Global.h"
#include "StepperBS.h"
#include "cbmpi.h"

/******************************************************************************
 * StepperBS updates x, y, dydx and gives interpolation values for dense output 
 * points.
 *
 * Input:
 *  x, y, dydx: changing variables in ODE.
 *  x0: array for dense output points.
 *  it1: last dense output point.
 *  densFlag: the flag for whether to do dense output.
 * Output:
 *  x, y, dydx, it1 are updated.
 *  hnext: the predicted next integration step.
 *  ydense: dense outputs at dense output points x0.
 ******************************************************************************/
void StepperBS(double *y, double *dydx, double *x, double *hnext, 
        double *x0, double *ydense, int *it1, bool densFlag)
{
    int nseq[IMAXX], cost[IMAXX]; /* Stepsize sequence and total steps */
    double coeff[IMAXX][IMAXX]; /* Coefficients used in extrapolation tableau */
    double errfac[2*IMAXX+2]; /* Used to compute dense interpolation error */
    double ysave[IMAXX][n]; /* ysave and fsave store values and derivatives */
    double fsave[IMAXX*(2*IMAXX+1)][n]; /* for dense output */
    double table[KMAXX][n];     /* Extrapolation tableau */
    double dens[(2*IMAXX+5)*n]; /* Stores quantities for dense interpolating 
                                   polynomial */
    int ipoint[IMAXX+1];        /* Keeps track of where values are stored 
                                   in fsave */
    double dydxnew[n];
    double h, htry, hdid;
    //double htry = 1.0E-5;
    double xold, mu;
    int i, j;
    int k_targ;

    htry = *hnext;
    xold = *x;
    mu = 0.0;
    for (i=0; i<IMAXX*(2*IMAXX+1); i++) {
        for (j=0; j<n; j++) {
            fsave[i][j] = 0.0;
        }
    }
    StepperBase(nseq, cost, hnext, &k_targ, coeff, errfac, ipoint, densFlag);
    *hnext = htry;
    //for (i=0; i<n; i++) {
    //    printf("i, dydx %d %20.14e\n", i, dydx[i]);
    //}
    step(htry, hnext, x, y, dydx, k_targ, nseq, table, cost, coeff,
         ipoint, errfac, dens, densFlag, dydxnew, ysave, fsave, &hdid);
    h = hdid;
  
    j = *it1;
    //printf("j                  %d\n", j);
    /* x0 is the array of dense output points. Find out all the dense
     * output points between xold and xold_h.
     */
    while ((x0[j]>=xold) && (x0[j]<=xold+h)) {
        for (i=0; i<n; i++) {
            ydense[(j-1)*nvar+i]=dense_out(i, x0[j], h, xold, mu, dens);
        }
        j++;
    }
    *it1 = j;
}

/******************************************************************************
 * StepperBase is used to initialize some of the variables in the code.
 ******************************************************************************/
void StepperBase(int *nseq, int *cost, double *hnext, int *k_targ, 
                 double coeff[][IMAXX], double *errfac, int *ipoint, bool densFlag)
{
    double logfact, ratio, e;
    int i, j, k, l, ip5, njadd;
    
    /* Stepsize sequence and total steps */
    if(densFlag) {
        for(i = 0; i < IMAXX; i++) nseq[i] = 4*i+2;
    }
    else {
        for(i = 0; i < IMAXX; i++) nseq[i] = 2*(i+1);
    }
    cost[0] = nseq[0] + 1;
    for(k = 0; k < KMAXX; k++) cost[k+1] = cost[k] + nseq[k+1];

    *hnext = -1.0E99;
    logfact = -log10(MAX(1.0e-12,rtol))*0.6+0.5;
    /* Initial estimate of optimal k */
    *k_targ = MAX(1,MIN(KMAXX-1,(int)logfact));

    /* Coefficients used in extrapolation tableau */
    for(k = 0; k < IMAXX; k++) {
        for(l = 0; l < k; l++) {
            ratio = (double)(nseq[k])/nseq[l];
            coeff[k][l] = 1.0/(ratio*ratio-1.0);
        }
    }

    /* Used to compute dense interpolation error */
    for(i = 0; i < 2*IMAXX+1; i++) {
        ip5 = i+5;
        errfac[i] = 1.0/(ip5*ip5);
        e = 0.5*sqrt(((double)(i+1))/ip5);
        for(j = 0; j <= i; j++) {
            errfac[i] *= e/(j+1);
        }
    }

    /* Keeps track of where values are stored in fsave */
    ipoint[0] = 0;
    for(i = 1; i <= IMAXX; i++) {
        njadd = 4*i - 2;
        if(nseq[i-1] > njadd) njadd++;
        ipoint[i] = ipoint[i-1] + njadd;
    }
}

/******************************************************************************
 * Attempts a step with stepsize htry. On output, y and x are replaced by their 
 * new values, hmid is the stepsize that was actually accomplished, and hnext 
 * is the estimated next stepsize. This is the actual stepper. It attempts a 
 * step, invokes the controller to decide whether to accept the step or try 
 * again with a smaller stepsize, and sets up the coefficients in case
 * dense output is needed between x and x+h.
 ******************************************************************************/
void step(double htry, double *hnext0, double *xpre, double *y, double *dydx, 
        int k_targ, int *nseq, double table[][n], int *cost, 
        double coeff[][IMAXX], int *ipoint, double *errfac, double *dens, 
        bool densFlag, double *dydxnew, double ysave[][n], double fsave[][n], 
        double *hdid)
{
    const double STEPFAC1=0.65, STEPFAC2=0.94, STEPFAC3=0.02, STEPFAC4=4.0,
                 KFAC1=0.8, KFAC2=0.9;
    static bool first_step=true, last_step=false;
    static bool forward, reject=false, prev_reject=false;
    double fac, h, hnew, hopt_int, err;
    double hopt[IMAXX], work[IMAXX];
    double ysav[n], yseq[n];
    double scale[n];
    double x, hnext;
    bool firstk;
    int i, k=0;

    x = *xpre;
    hnext = *hnext0;
    work[0] = 0.0;
    h = htry;
    forward = h>0 ? true : false;
    for(i = 0; i < n; i++) ysav[i] = y[i]; /* Saving the starting values. */
    /* h gets reset in Odeint for the last step. */
    if(h != hnext && !first_step) {
        last_step=true;
    }
    if(reject) {
        prev_reject = true;
        last_step = false;
    }
    reject = false;
    firstk = true;
    hnew = fabs(h);

    /* Restart here if interpolation error is too big. Loop until step accepted. */
interp_error:
    while(firstk || reject) {
        h = forward ? hnew : -hnew;
        firstk = false;
        reject = false;
        if(fabs(h) <= fabs(x)*EPS) {
            printf("step size underflow in StepperBS\n");
        }
        // Initialize counter for saving stuff.
        int ipt = -1;
        // Evaluate the sequence of modified midpoint integrations.
        for(k = 0; k <= k_targ+1; k++) {
            //printf("k = %d\n", k);
            updatey(x, ysav, dydx, h, k, ipt, nseq, ysave, fsave, yseq, densFlag);
            if(k == 0) {
                for (i=0; i<n; i++) {
                    y[i] = yseq[i];
                }
                //y = yseq;
            }
            else {
                for (i=0; i< n; i++) table[k-1][i] = yseq[i];
            }
            if(k != 0) {
                //printf("1 %20.14e\n", y[0]);
                polyextr(k, table, y, coeff);
                //printf("2 %20.14e\n", y[0]);
                err = 0.0;
                // Compute normalized error estimate err.
                for(i = 0; i < n; i++) {
                    scale[i]=atol+rtol*MAX(fabs(ysav[i]), fabs(y[i]));
                    err += pow((y[i]-table[0][i])/scale[i], 2.0);
                }
                err = sqrt(err/n);
                // Compute optimal stepsize for this order.
                double expo = 1.0/(2*k+1); // Compute optimal stepsize for this order.
                double facmin = pow(STEPFAC3,expo);
                if(err == 0.0)
                    fac = 1.0/facmin;
                else {
                    fac = STEPFAC2/pow(err/STEPFAC1,expo);
                    fac = MAX(facmin/STEPFAC4,MIN(1.0/facmin,fac));
                }
                hopt[k] = fabs(h*fac);
                work[k] = cost[k]/hopt[k];
                if((first_step || last_step) && err <= 1.0) {
                    //printf("1st 2nd %d %d\n", first_step, last_step);
                    //printf("error2 %20.14e\n", err);
                    break;
                }
                if(k == k_targ-1 && !prev_reject && !first_step && !last_step) {
                    if(err <= 1.0) // Converged within order window.
                        break;
                    else if(err > pow(nseq[k_targ]*nseq[k_targ+1]/(nseq[0]*nseq[0]), 2.0)) {
                        reject = true;
                        k_targ = k;
                        if (k_targ > 1 && work[k-1] < KFAC1*work[k]) k_targ--;
                        hnew = hopt[k_targ];
                        break;
                    }
                }
                if(k == k_targ) {
                    if(err <= 1.0)
                        break;
                    else if(err > pow((double)nseq[k+1]/nseq[0], 2.0)) {
                        reject = true;
                        if (k_targ > 1 && work[k-1] < KFAC1*work[k]) k_targ--;
                        hnew = hopt[k_targ];
                        break;
                    }
                }
                if(k == k_targ + 1) {
                    if(err > 1.0) {
                        reject = true;
                        if(k_targ > 1 && work[k_targ-1] < KFAC1*work[k_targ]) k_targ--;
                        hnew = hopt[k_targ];
                    }
                    break;
                }
            }
        }
        // Go back and try next k.
        // Arrive here from any break in for loop.
        if(reject) prev_reject = true;
    }

    //printf("hnext: %20.14e\n", hnew);

    // Go back if step was rejected.
    // Used for start of next step and in dense output.
    derivs(x+h, y, dydxnew);
    if (densFlag) {
        prepare_dense(h, dydxnew, ysav, scale, k, &err, ipoint, 
                      errfac, y, dydx, nseq, ysave, fsave, dens);
        hopt_int = h/MAX(pow(err, 1.0/(2*k+3)), 0.01);
        //printf("error: %20.14e\n", err);
        if(err > 10.0) {
            hnew = fabs(hopt_int);
            reject = true;
            prev_reject = true;
            goto interp_error;
        }
    }

    for (i=0; i<n; i++) {
        dydx[i] = dydxnew[i];
    }
    //dydx = dydxnew;
    //xold = x;
    x += h;
    *xpre = x;
    *hdid = h;
    //printf("hdid %20.14e\n", *hdid);
    first_step = false;
    int kopt;
    if(k == 1)
        kopt = 2;
    else if(k <= k_targ) {
        kopt = k;
        if(work[k-1] < KFAC1*work[k])
            kopt = k-1;
        else if(work[k] < KFAC2*work[k-1]) {
            kopt = MIN(k+1,KMAXX-1);
        }
    } else {
        kopt = k-1;
        if(k > 2 && work[k-2] < KFAC1*work[k-1])
            kopt = k-2;
        if(work[k] < KFAC2*work[kopt])
            kopt = MIN(k, KMAXX-1);
    }
    if(prev_reject) {
        k_targ = MIN(kopt,k);
        hnew = MIN(fabs(h), hopt[k_targ]);
        prev_reject = false;
    }
    else {
        if(kopt <= k)
            hnew = hopt[kopt];
        else {
            if(k < k_targ && work[k] < KFAC2*work[k-1])
                hnew = hopt[k]*cost[kopt+1]/cost[k];
            else
                hnew = hopt[k]*cost[kopt]/cost[k];
        }
        k_targ = kopt;
    }
    if(densFlag)
        hnew = MIN(hnew,fabs(hopt_int));
    if(forward) {
        hnext = hnew;
    }
    else {
        hnext = -hnew;
    }

    *hnext0 = hnext;
    //printf("y final %20.14e\n", y[0]);
}

/***********************************************************************************************
 * Given values for n variables y[0..n-1] and their derivatives dydx[0..n-1] known at x, use the
 * Bulirsch-Stoer method to advance the solution over an interval h.
 * z0 = f(x)
 * z1 = z0+h*f(x,z0)
 * z(m+1) = z(m-1) + 2*h*f(x+m*h,z(m))
 * y(x+H) ~= yn = 0.5*(z(n)+z(n-1)+h*f(x+H,z(n)))
 ***********************************************************************************************/
void updatey(double x, double *y, double *dydx, const double htot, const int k, int ipt, 
             int *nseq, double ysave[][n], double fsave[][n], double *yend, bool densFlag)
{
    double ym[n], yn[n];
    double xnew, h2, swap;
    double h;
    int nstep = nseq[k];
    int i, nn;

    h = htot/nstep;
    for(i = 0; i < n; i++) {
        ym[i] = y[i];
        yn[i] = y[i] + h*dydx[i];
    }
    //printf("------------------ h ------------ %20.14e\n", h);
    xnew = x + h;
    // Use yend for temporary storage of derivatives.
    derivs(xnew, yn, yend);
    h2 = 2.0*h;
    // for dense output
    for(nn = 1; nn < nstep; nn++) {
        if(densFlag && nn == nstep/2) {
            for(i = 0; i < n; i++) {
                ysave[k][i] = yn[i];
            }
        }
        if(densFlag && fabs(nn-nstep/2) <= 2*k+1) {
            ipt++;
            for(i = 0; i < n; i++)
                fsave[ipt][i] = yend[i];
        }
        for(i = 0; i < n; i++) {
            swap = ym[i] + h2*yend[i];
            ym[i] = yn[i];
            yn[i] = swap;
        }
        xnew += h;
        derivs(xnew, yn, yend);
    }
    if(densFlag && nstep/2 <= 2*k+1) {
        ipt++;
        for (i = 0; i < n; i++)
            fsave[ipt][i] = yend[i];
    }
    
    // final point of one step
    for(i = 0; i < n; i++) {
        yend[i] = 0.5*(ym[i]+yn[i]+h*yend[i]);
        //printf("y1, y2 = %20.14e %20.14e\n", yend[i], y[i]);
    }
}

/***********************************************************************************************
 * Polynomial extrapolation
 ***********************************************************************************************/
void polyextr(const int k, double table[][n], double *last, double coeff[][IMAXX])
{
    int i, j, l;
    l = n;
    for (j = k-1; j > 0; j--) {
        for(i = 0; i < l; i++) {
            table[j-1][i] = table[j][i] + coeff[k][j]*(table[j][i]-table[j-1][i]);
        }
    }
    //printf("coeff %20.14e\n", coeff[k][0]);
    for (i = 0; i < l; i++) {
        //printf("table %d %20.14e %20.14e\n", i, table[0][i], last[i]);
        //printf("last%d %20.14e\n", i, last[i]);
        last[i] = table[0][i] + coeff[k][0]*(table[0][i]-last[i]);
        //printf("last%d %20.14e\n", i, last[i]);
    }
}

/***********************************************************************************************
 * Derivatives used in the calculation.  For both of x and v.
 * y is array with elements x, y, z, vx, vy, vz. 
 * yderiv is array with elements vx, vy, vz, ax, ay, az. 
 * All v here is actually gamma = v/c, so a is d(gamma) / dt.
 ***********************************************************************************************/
void derivs(double t, double *y, double *yderiv)
{
    struct emfields emf;
    double Cx, Cy, Cz;
    const double cl=c0/L0;
    double beta2, gama, igama3;
    int i;
    // for x
    for(i = 0; i < n/2; i++) {
        yderiv[i] = y[i+n/2]*cl;
        //printf("beta %20.14e\n", yderiv[i]);
    }
    // for v, actually beta
    emf.Bx = 0.0; emf.By = 0.0; emf.Bz = 0.0;
    emf.Ex = 0.0; emf.Ey = 0.0; emf.Ez = 0.0;
    Cx = 0.0; Cy = 0.0; Cz = 0.0;
    get_emf(y[0], y[1], y[2], t, &emf);
	cross_product(y[3], y[4], y[5], emf.Bx, emf.By, emf.Bz, &Cx, &Cy, &Cz);
    beta2 = y[3]*y[3] + y[4]*y[4] + y[5]*y[5];
    //printf("=========beta2= %20.14e\n", beta2);
    gama = 1.0/sqrt(1.0-beta2);
    igama3 = 1.0/(gama*gama*gama);
    yderiv[3] = igama3*charge_mass*(emf.Ex+Cx);
    yderiv[4] = igama3*charge_mass*(emf.Ey+Cy);
    yderiv[5] = igama3*charge_mass*(emf.Ez+Cz);
    //for(i = n/2; i < n; i++) {
    //    printf("beta %20.14e\n", yderiv[i]);
    //}
}

/***********************************************************************************************
 * dense output interpretation. Compute coefficients of the dense interpolation formula. On
 * input, y[0...n*(imit+5)-1] contains the dens array from prepare-dense. On output these
 * coefficients have been updated to the required values.
 ***********************************************************************************************/
void dense_interp(double *y, const int imit)
{
    double y0, y1, yp0, yp1, ydiff, aspl, bspl, ph0, ph1, ph2, ph3, fac1, fac2;
    double a[31];
    int i, im;
    for (i = 0; i < n; i++) {
        y0 = y[i];
        y1 = y[2*n+i];
        yp0 = y[n+i];
        yp1 = y[3*n+i];
        ydiff = y1 - y0;
        aspl = -yp1 + ydiff;
        bspl = yp0 - ydiff;
        y[n+i] = ydiff;
        y[2*n+i] = aspl;
        y[3*n+i] = bspl;
        if (imit < 0) continue;
        ph0 = (y0+y1)*0.5 + 0.125*(aspl+bspl);
        ph1 = ydiff + (aspl-bspl)*0.25;
        ph2 = -(yp0-yp1);
        ph3 = 6.0 * (bspl-aspl);
        if (imit >= 1) {
            a[1] = 16.0*(y[5*n+i]-ph1);
            if (imit >= 3) {
                a[3] = 16.0*(y[7*n+i]-ph3+3*a[1]);
                for (im = 5; im <= imit; im += 2) {
                    fac1 = im*(im-1) / 2.0;
                    fac2 = fac1*(im-2)*(im-3)*2.0;
                    a[im] = 16.0*(y[(im+4)*n+i]+fac1*a[im-2]-fac2*a[im-4]);
                }
            }
        }
        a[0] = (y[n*4+i]-ph0)*16.0;
        if (imit >= 2) {
            a[2] = (y[n*6+i]-ph2+a[0])*16.0;
            for (im = 4; im <= imit; im += 2) {
                fac1 = im*(im-1) / 2.0;
                fac2 = im*(im-1)*(im-2)*(im-3);
                a[im] = 16.0*(y[(im+4)*n+i]+fac1*a[im-2]-fac2*a[im-4]);
            }
        }
        for (im = 0; im <= imit; im++) {
            y[n*(im+4)+i] = a[im];
        }
    }
}

/***********************************************************************************************
 * dense output preparation to set up the dense output quantities.
 * Store coefficients of interpolating polynomial for dense output in dens array. Input stepsize
 * h, derivative at end of interval dydexnew[0...n-1], function at beginning of interval
 * ysav[0...n-1], scale factor atol+|y|*rtol in scale[0...n-1], and column k in which convergence
 * was achieved. Output interpolation error in error.
 ***********************************************************************************************/
void prepare_dense(const double h, double *dydxnew, double *ysav, double *scale, const int k, 
                   double *error, int *ipoint, double *errfac, double *y, double *dydx, 
                   int *nseq, double ysave[][n], double fsave[][n], double *dens)
{
    double dblenj, factor, facnj;
    int i, j, l, mu, kmi, kk, kbeg, ipt, lbeg, lend;
    mu = 2*k-1;
    for (i = 0; i < n; i++) {
        //printf("============dydx, dydxnew= %lf %lf\n", dydx[i], dydxnew[i]);
        dens[i] = ysav[i];
        dens[n+i] = h*dydx[i];
        dens[2*n+i] = y[i];
        dens[3*n+i] = dydxnew[i]*h;
    }
    for (j = 1; j <= k; j++) {
        dblenj = nseq[j];
        for (l = j; l >= 1; l--) {
            factor = pow(dblenj/nseq[l-1],2.0) - 1.0;
            for (i = 0; i < n; i++) {
                ysave[l-1][i] = ysave[l][i] + (ysave[l][i]-ysave[l-1][i])/factor;
            }
        }
    }
    for (i = 0; i < n; i++) {
        dens[4*n+i] = ysave[0][i];
    }
    for (kmi = 1; kmi <= mu; kmi++) {
        kbeg = (kmi-1) / 2;
        for (kk = kbeg; kk <= k; kk++) {
            facnj = pow(nseq[kk]/2.0, kmi-1);
            ipt = ipoint[kk+1] - 2*kk + kmi - 3;
            for (i = 0; i < n; i++) {
                ysave[kk][i] = fsave[ipt][i]*facnj;
                //printf("============fsave= %20.14e \n", fsave[ipt][i]);
            }
        }
        for (j = kbeg+1; j <= k; j++) {
            dblenj = nseq[j];
            for (l = j; l >= kbeg+1; l--) {
                factor = pow(dblenj/nseq[l-1],2.0) - 1.0;
                for (i = 0; i < n; i++) {
                    ysave[l-1][i] = ysave[l][i] + (ysave[l][i]-ysave[l-1][i])/factor;
                }
            }
        }
        for (i = 0; i < n; i++) {
            dens[(kmi+4)*n+i] = ysave[kbeg][i] * h;
        }
        if (kmi == mu) continue;
        // Compute differences.
        for (kk = kmi/2; kk <= k; kk++) {
            lbeg = ipoint[kk+1] - 1;
            lend = ipoint[kk] + kmi;
            if (kmi == 1) lend += 2;
            for (l = lbeg; l >= lend; l -= 2) {
                for (i = 0; i < n; i++) {
                    fsave[l][i] = fsave[l][i] - fsave[l-2][i];
                }
            }
            if (kmi == 1) {
                l = lend - 2;
                for (i = 0; i < n; i++) {
                    fsave[l][i] = fsave[l][i] - dydx[i];
                }
            }
        }
        for (kk = kmi/2; kk <= k; kk++) {
            lbeg = ipoint[kk+1] - 2;
            lend = ipoint[kk] + kmi + 1;
            for (l = lbeg; l >= lend; l -= 2) {
                for (i = 0; i < n; i++) {
                    fsave[l][i] = fsave[l][i] - fsave[l-2][i];
                }
            }
        }
    }
    // Compute the interpolation coefficients in dens. 
    // Estimate the interpolation error.
    dense_interp(dens, mu);
    *error = 0.0;
    if (mu >= 1) {
        for (i = 0; i < n; i++) {
            *error += pow(dens[(mu+4)*n+i]/scale[i], 2);
            //printf("dens %20.14e %20.14e\n", dens[(mu+4)*n+i], scale[i]);
        }
        *error = sqrt(*error/n)*errfac[mu-1];
    }
}

/******************************************************************************
 * dense output, using the coefficients stored by the previous routine to
 * evaluate the solution at an arbitrary point. Evaluate interpolating polynomial 
 * for y[i] at location x, where xold <= x <= xold+h.
 ******************************************************************************/
double dense_out(const int i, const double x, const double h, double xold, int mu, double *dens)
{
    double theta, theta1, theta05, t4, c;
    double yinterp;
    int j;
    theta = (x-xold) / h;
    theta1 = 1.0 - theta;
    yinterp = dens[i] + theta*(dens[n+i]+theta1*(dens[2*n+i]*theta+dens[3*n+i]*theta1));
    if (mu < 0)
        return yinterp;
    theta05 = theta - 0.5;
    t4 = pow(theta*theta1,2);
    c = dens[n*(mu+4)+i];
    for (j = mu; j > 0; j--) {
        c = dens[n*(j+3)+i] + c*theta05/j;
    }
    yinterp += t4 * c;
    return yinterp;
}
