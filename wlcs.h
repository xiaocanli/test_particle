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
void getWireB(double x, double y, double z, int iconf, struct bfields *bmf);
void getWire_Bessel(double x, double y, double z, int iconf, double t, 
        struct emfields *emf);
void getLoopB(double x, double y, double z, int iconf, struct bfields *bmf);
void ellint(double k0, double *elle, double *ellf);
void transform(double *Ax, double *Ay, double *Az, 
        double cos_euler1, double sin_euler1, 
        double cos_euler2, double sin_euler2);
void reverse_tran(double *Ax, double *Ay, double *Az, 
        double cos_euler1, double sin_euler1, 
        double cos_euler2, double sin_euler2);
void read_wlcs();
void getemf_wlcs(double x, double y, double z, double t, 
        struct emfields *emf_tot);
