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
#define KMAXX 8
#define IMAXX (KMAXX+1)
#define EPS 2.22045e-16
#define atol 1.0E-5
#define rtol 1.0E-14
#define n 6

void StepperBase(int *nseq, int *cost, double *hnext, int *k_targ, 
        double coeff[][IMAXX], double *errfac, int *ipoint, bool densFlag);
void step(double htry, double *hnext0, double *xpre, double *y, 
        double *dydx, int k_targ, int *nseq,
        double table[][n], int *cost, double coeff[][IMAXX],
        int *ipoint, double *errfac, double *dens, bool densFlag, 
        double *dydxnew, double ysave[][n], double fsave[][n], double *hdid);
void updatey(double x, double *y, double *dydx, const double htot, 
        const int k, int ipt, int *nseq, double ysave[][n], 
        double fsave[][n], double *yend, bool densFlag);
void polyextr(const int k, double table[][n], double *last, 
        double coeff[][IMAXX]);
void prepare_dense(const double h, double *dydxnew, double *ysav, 
        double *scale, const int k, double *error, int *ipoint, 
        double *errfac, double *y, double *dydx, int *nseq, 
        double ysave[][n], double fsave[][n], double *dens);
double dense_out(const int i, const double x, const double h, 
        double xold, int mu, double *dens);
void dense_interp(double *y, const int imit);
void StepperBS(double *y, double *dydx, double *x, double *hnext, 
        double *x0, double *ydense, int *it1, bool densFlag);
void derivs(double t, double *y, double *yderiv);
