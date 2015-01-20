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
#include "Global.h"
#include "wlcs.h"

struct wlcs *config;
int nwlcs, iBessel;

/******************************************************************************
 * Read the configuration data of a wire-loop current system (WLCS).
 * Sin and Cos functions are calculated here since they  will be re-used. 
 *
 * Output:
 *  config: configuration array for wire-loop current system.
 ******************************************************************************/
void read_wlcs()
{
	FILE *fp;
    double theta, phi, alpha, beta, omega;
    char buff[200];
    int i, msg;
    fp = fopen("init.dat", "r");
    if (fgets(buff, 200, fp) != NULL) {
        puts(buff);
    } else {
        printf("Error: hit the end the of file.\n");
        exit(1);
    }
    nwlcs = 1;
    msg = fscanf(fp, "Number of WLCSs:%d\n", &nwlcs);
    if (msg == 1) {
        printf("Number of WLCSs: %d\n", nwlcs);
    } else {
        printf("Failed to read the number of WLCSs.\n");
        exit(1);
    }
    fclose(fp);

    config = (struct wlcs*)malloc(sizeof(struct wlcs)*nwlcs);

	fp = fopen("init8_3l.dat","r");
    if (fgets(buff, 200, fp)) {
        
    };
	for(i = 0; i<nwlcs; i++) {
        msg = fscanf(fp, "%lf", &(config[i].cur_wr));
        if (msg != 1) {
            printf("Failed to read wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].x_wr));
        if (msg != 1) {
            printf("Failed to read x position of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].y_wr));
        if (msg != 1) {
            printf("Failed to read the y position of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].z_wr));
        if (msg != 1) {
            printf("Failed to read the z position of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &theta);
        if (msg != 1) {
            printf("Failed to read theta Euler angle of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &phi);
        if (msg != 1) {
            printf("Failed to read phi Euler angle of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].t0));
        if (msg != 1) {
            printf("Failed to read t0 of wire current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].cur_lp));
        if (msg != 1) {
            printf("Failed to read the loop current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].r_lp));
        if (msg != 1) {
            printf("Failed to read the radius of the loop current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].x_lp));
        if (msg != 1) {
            printf("Failed to read the x position of the loop center.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].y_lp));
        if (msg != 1) {
            printf("Failed to read the y position of the loop center.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &(config[i].z_lp));
        if (msg != 1) {
            printf("Failed to read the z position of the loop center.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &alpha);
        if (msg != 1) {
            printf("Failed to read alpha Euler angle of the loop current.\n");
            exit(1);
        }
		msg = fscanf(fp, "%lf", &beta);
        if (msg != 1) {
            printf("Failed to read beta Euler angle of the loop current.\n");
            exit(1);
        }

        /* avoiding re-calculation for magnetic field */
		config[i].costheta_wr = cos(theta);
		config[i].sintheta_wr = sin(theta);
		config[i].cosphi_wr = cos(phi);
		config[i].sinphi_wr = sin(phi);
		config[i].cosalpha_lp = cos(alpha);
		config[i].sinalpha_lp = sin(alpha);
		config[i].cosbeta_lp = cos(beta);
		config[i].sinbeta_lp = sin(beta);
	}

	msg = fscanf(fp, "%d", &iBessel);
    if (msg != 1) {
        printf("Failed to read iBessel.\n");
        exit(1);
    }
	msg = fscanf(fp, "%lf", &(omega));
    if (msg != 1) {
        printf("Failed to read omega (frequency of the wire currents).\n");
        exit(1);
    }

    for(i = 0; i < nwlcs; i++) {
        config[i].omega = omega;
    }

//    for(i = 0; i < nwlcs; i++){
//		printf("%lf ", config[i].cur_wr);
//		printf("%lf ", config[i].x_wr);
//		printf("%lf ", config[i].y_wr);
//		printf("%lf ", config[i].z_wr);
//		printf("%lf ", config[i].sintheta_wr);
//		printf("%lf ", config[i].costheta_wr);
//		printf("%lf ", config[i].sinphi_wr);
//		printf("%lf ", config[i].cosphi_wr);
//		printf("%lf ", config[i].t0);
//		printf("%lf ", config[i].cur_lp);
//		printf("%lf ", config[i].r_lp);
//		printf("%lf ", config[i].x_lp);
//		printf("%lf ", config[i].y_lp);
//		printf("%lf ", config[i].z_lp);
//		printf("%lf ", config[i].sinalpha_lp);
//		printf("%lf ", config[i].cosalpha_lp);
//		printf("%lf ", config[i].sinbeta_lp);
//		printf("%lf ", config[i].cosbeta_lp);
//        printf("%lf ", config[i].omega);
//        printf("\n");
//    }
	fclose(fp);
}

/******************************************************************************
 * Calculate the magnetic field of one straight wire in time-dependent condition.
 *
 * Input:
 *  x, y, z, t: spatial positions and time
 *  iconf: the id of wlcs
 * Output:
 *  emf: electromagnetic fields.
 ******************************************************************************/
void getWire_Bessel(double x, double y, double z, 
        int iconf, double t, struct emfields *emf)
{
	double xp, yp, zp;
	double x0, rho, const0;
	double bj0, bj1, by0, by1; // Bessel functions
	double wireB, wireE;
	double phi1, curI, omega;
    double t1;

//    printf("x, y, z: %lf %lf %lf %lf\n", x, y, z, t);

    emf->Bx = 0.0; emf->By = 0.0; emf->Bz = 0.0;
    emf->Ex = 0.0; emf->Ey = 0.0; emf->Ez = 0.0;

	xp = x - config[iconf].x_wr;
	yp = y - config[iconf].y_wr;
	zp = z - config[iconf].z_wr;
    t1 = config[iconf].t0 + t;
    curI = config[iconf].cur_wr;
    omega = config[iconf].omega;

    transform(&xp, &yp, &zp, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr,
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
	rho = sqrt(xp*xp + yp*yp); 
    x0 = 2.32E-2 * rho * omega; // page 9
	//bj0 = bessj0(x0);
	//by0 = bessy0(x0);
	//bj1 = bessj1(x0);
	//by1 = bessy1(x0);
    bj0 = j0(x0);
	by0 = y0(x0);
    bj1 = j1(x0);
	by1 = y1(x0);

	const0 = (curI*I0) / (5.0*rho*L0) * M_PI/2.0;

//---------------------------------- Bx, By, Bz ----------------------------------
    wireB = const0 * x0 * (-by1*cos(omega * t1) + bj1 * sin(omega * t1));

    if(rho == 0.0){
		phi1 = M_PI/2.0;
	}
	else if(yp >= 0.0) {
		phi1 = acos(xp/rho) + M_PI/2.0;
	}
	else {
		phi1 = 2.0*M_PI - acos(xp/rho) + M_PI/2.0;
	}
		
	emf->Bx = wireB * cos(phi1);
	emf->By = wireB * sin(phi1);
	emf->Bz = 0.0;
		
    reverse_tran(&emf->Bx, &emf->By, &emf->Bz, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr, 
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);

//---------------------------------- Ex, Ey, Ez ----------------------------------
	wireE = -const0 * x0 * (by0*sin(omega * t1) + bj0*cos(omega * t1));

	emf->Ex = 0.0;
	emf->Ey = 0.0;
	emf->Ez = wireE;

    reverse_tran(&emf->Ex, &emf->Ey, &emf->Ez, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr,
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
}

/******************************************************************************
 * Calculate the magnetic field of one straight wire in time-independent condition.
 ******************************************************************************/
void getWireB(double x, double y, double z, int iconf, struct bfields *bmf)
{
	double xp, yp, zp;
	double rho, cos_phi0, phi0;
	double Bphi, curI;

    bmf->Bx = 0.0; bmf->By = 0.0; bmf->Bz = 0.0; // initilization
	
    xp = x - config[iconf].x_wr;
	yp = y - config[iconf].y_wr;
	zp = z - config[iconf].z_wr;
    curI = config[iconf].cur_wr;

    transform(&xp, &yp, &zp, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr, 
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
	rho = sqrt(xp*xp + yp*yp); // distance

	if(rho == 0.0) {
		printf("%s\n","=== trying to calculate B at the wire ===");
		return;
	}
	else {
		cos_phi0 = xp / rho;
		if(yp >= 0.0)
			phi0 = acos(cos_phi0);
		else
			phi0 = 2.0*M_PI - acos(cos_phi0);
	}

	Bphi = 1.0/5.0 * (curI*I0) / (rho*L0); // phi component

	bmf->Bx = Bphi * cos(phi0 + M_PI/2.0);
	bmf->By = Bphi * sin(phi0 + M_PI/2.0);
	bmf->Bz = 0.0;

    reverse_tran(&bmf->Bx, &bmf->By, &bmf->Bz, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr, 
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
}

/******************************************************************************
 * Calculate the magnetic field of one loop in time-independent condition.
 ******************************************************************************/
void getLoopB(double x, double y, double z, int iconf, struct bfields *bmf)
{
	double xp, yp, zp, curI, al;
	double k_tmp, E_tmp, F_tmp;
	double x_tmp, rho_tmp, k2_tmp, BpRho4, BpRho3;
    double tmp1, tmp2, tmp3, tmp4;
    double rhop1, rhop3;
    double cosalpha, sinalpha, cosbeta, sinbeta;
	
    bmf->Bx = 0.0; bmf->By = 0.0; bmf->Bz = 0.0; // initilization

    // transformation from different coordinate systems	
	xp = (x - config[iconf].x_lp) * L0;
	yp = (y - config[iconf].y_lp) * L0;
	zp = (z - config[iconf].z_lp) * L0;
    al = config[iconf].r_lp * L0;
    curI = config[iconf].cur_lp * I0;
    cosalpha = config[iconf].cosalpha_lp;
    sinalpha = config[iconf].sinalpha_lp;
    cosbeta = config[iconf].cosbeta_lp;
    sinbeta = config[iconf].sinbeta_lp;

    transform(&xp, &yp, &zp, cosbeta, sinbeta, cosalpha, sinalpha);

	rhop3 = sqrt(xp*xp + yp*yp);
    rhop1 = sqrt((rhop3 + al)*(rhop3 + al) + zp * zp);

    // The magnetic field (Bpx3, Bpy3, Bpz3) is obtained in
    // the frame defined by (xp, yp, zp)
	k_tmp = sqrt(4.0*al*rhop3) / rhop1;
    ellint(k_tmp, &E_tmp, &F_tmp);
	x_tmp = 2.0 * al / sqrt(al*al + zp*zp);
	rho_tmp = rhop3 / al;
	k2_tmp = x_tmp * sqrt(rho_tmp) / sqrt(1.0 + 0.25 * x_tmp * x_tmp * 
            (rho_tmp * rho_tmp + 2.0 * rho_tmp));
	BpRho4 = 0.2 * curI * zp / (2.0 * al * al) * 
            (3.0*M_PI/32.0) * pow(k2_tmp, 5) * 
            pow(1.0/rho_tmp, 3.0/2.0) / (1.0 - k2_tmp*k2_tmp);

    tmp1 = 0.2 * curI / rhop1;
    tmp2 = rhop3 * rhop3 + zp * zp;
    tmp3 = al * al;
    tmp4 = E_tmp / (tmp2 + tmp3 - 2.0 * rhop3 * al);

	BpRho3 = tmp1 * zp * (-F_tmp + (tmp2 + tmp3) * tmp4) / rhop3;

	if(rhop3 == 0.0){
		BpRho3 = BpRho4;
	}

    BpRho3 /= rhop3; // for saving computing time
	bmf->Bx = BpRho3 * xp;
	bmf->By = BpRho3 * yp;
	bmf->Bz = tmp1 * (F_tmp - (tmp2 - tmp3) * tmp4);
        
    reverse_tran(&bmf->Bx, &bmf->By, &bmf->Bz, 
            cosbeta, sinbeta, cosalpha, sinalpha);
}

/******************************************************************************
 * Elliptic integral using Toshio Fukushima's method.
 * Fast computation of complete elliptic integrals and Jacobian elliptic functions,
 * Toshio Fukushima, Celest Mech Dyn Astr (2009) 105:305–328
 * DOI 10.1007/s10569-009-9228-z
 *
 * Author: Xiaocan Li
 * Date: Oct. 22. 2012
 ******************************************************************************/
void ellint(double k, double *elle, double *ellf)
{
    double m, Kj[20], Ej[20], qj[15], Hj[11];
    double ellf1 = 0.0; // associate integrals
    double q1 = 0.0, Hm = 0.0; // used when m is [0.9, 1.0]
    int i, j;
    double powj = 1.0;
    m = k * k; // elliptic parameter
    *elle = 0.0; *ellf = 0.0;

    if(m < 0.8 ) {
        i = (int)(m/0.1);
    }
    else if(m >= 0.8 && m < 0.85) {
        i = 8;
    }
    else if(m >= 0.85 && m < 0.9) {
        i = 9;
    }
    else {
        i = 10;
    }
//    printf("k = %f\n", k);
//    printf("%d\n", i);
// Basically, we have to use Taylor series expansion for different
// intervals of m. The coefficients are given in Toshio's article.
    switch(i) {
        case 0:
            Kj[0] = 1.591003453790792180; Ej[0] = +1.550973351780472328;
            Kj[1] = 0.416000743991786912; Ej[1] = -0.400301020103198524;
            Kj[2] = 0.245791514264103415; Ej[2] = -0.078498619442941939;
            Kj[3] = 0.179481482914906162; Ej[3] = -0.034318853117591992;
            Kj[4] = 0.144556057087555150; Ej[4] = -0.019718043317365499;
            Kj[5] = 0.123200993312427711; Ej[5] = -0.013059507731993309;
            Kj[6] = 0.108938811574293531; Ej[6] = -0.009442372874146547;
            Kj[7] = 0.098853409871592910; Ej[7] = -0.007246728512402157;
            Kj[8] = 0.091439629201749751; Ej[8] = -0.005807424012956090;
            Kj[9] = 0.085842591595413900; Ej[9] = -0.004809187786009338;
            Kj[10] = 0.081541118718303215;
            for(j = 0; j <= 9; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.05;
            }
            *ellf += Kj[10] * powj;
            break;
        case 1:
            Kj[0] = 1.635256732264579992; Ej[0] = +1.510121832092819728;
            Kj[1] = 0.471190626148732291; Ej[1] = -0.417116333905867549;
            Kj[2] = 0.309728410831499587; Ej[2] = -0.090123820404774569;
            Kj[3] = 0.252208311773135699; Ej[3] = -0.043729944019084312;
            Kj[4] = 0.226725623219684650; Ej[4] = -0.027965493064761785;
            Kj[5] = 0.215774446729585976; Ej[5] = -0.020644781177568105;
            Kj[6] = 0.213108771877348910; Ej[6] = -0.016650786739707238;
            Kj[7] = 0.216029124605188282; Ej[7] = -0.014261960828842520;
            Kj[8] = 0.223255831633057896; Ej[8] = -0.012759847429264803;
            Kj[9] = 0.234180501294209925; Ej[9] = -0.011799303775587354;
            Kj[10] = 0.248557682972264071; Ej[10] = -0.011197445703074968;
            Kj[11] = 0.266363809892617521;
            for(j = 0; j <= 10; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.15;
            }
            *ellf += Kj[11] * powj;
            break;
        case 2:
            Kj[0] = 1.685750354812596043; Ej[0] = +1.467462209339427155;
            Kj[1] = 0.541731848613280329; Ej[1] = -0.436576290946337775;
            Kj[2] = 0.401524438390690257; Ej[2] = -0.105155557666942554;
            Kj[3] = 0.369642473420889090; Ej[3] = -0.057371843593241730;
            Kj[4] = 0.376060715354583645; Ej[4] = -0.041391627727340220;
            Kj[5] = 0.405235887085125919; Ej[5] = -0.034527728505280841;
            Kj[6] = 0.453294381753999079; Ej[6] = -0.031495443512532783;
            Kj[7] = 0.520518947651184205; Ej[7] = -0.030527000890325277;
            Kj[8] = 0.609426039204995055; Ej[8] = -0.030916984019238900;
            Kj[9] = 0.724263522282908870; Ej[9] = -0.032371395314758122;
            Kj[10] = 0.871013847709812357; Ej[10] = -0.034789960386404158;
            Kj[11] = 1.057652872753547036;
            for(j = 0; j <= 10; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.25;
            }
            *ellf += Kj[11] * powj;
            break;
        case 3:
            Kj[0] = 1.744350597225613243; Ej[0] = +1.422691133490879171;
            Kj[1] = 0.634864275371935304; Ej[1] = -0.459513519621048674;
            Kj[2] = 0.539842564164445538; Ej[2] = -0.125250539822061878;
            Kj[3] = 0.571892705193787391; Ej[3] = -0.078138545094409477;
            Kj[4] = 0.670295136265406100; Ej[4] = -0.064714278472050002;
            Kj[5] = 0.832586590010977199; Ej[5] = -0.062084339131730311;
            Kj[6] = 1.073857448247933265; Ej[6] = -0.065197032815572477;
            Kj[7] = 1.422091460675497751; Ej[7] = -0.072793895362578779;
            Kj[8] = 1.920387183402304829; Ej[8] = -0.084959075171781003;
            Kj[9] = 2.632552548331654201; Ej[9] = -0.102539850131045997;
            Kj[10] = 3.652109747319039160; Ej[10] = -0.127053585157696036;
            Kj[11] = 5.115867135558865806; Ej[11] = -0.160791120691274606;
            Kj[12] = 7.224080007363877411;
            for(j = 0; j <= 11; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.35;
            }
            *ellf += Kj[12] * powj;
            break;
        case 4:
            Kj[0] = 1.813883936816982644; Ej[0] = +1.375401971871116291;
            Kj[1] = 0.763163245700557246; Ej[1] = -0.487202183273184837;
            Kj[2] = 0.761928605321595831; Ej[2] = -0.153311701348540228;
            Kj[3] = 0.951074653668427927; Ej[3] = -0.111849444917027833;
            Kj[4] = 1.315180671703161215; Ej[4] = -0.108840952523135768;
            Kj[5] = 1.928560693477410941; Ej[5] = -0.122954223120269076;
            Kj[6] = 2.937509342531378755; Ej[6] = -0.152217163962035047;
            Kj[7] = 4.594894405442878062; Ej[7] = -0.200495323642697339;
            Kj[8] = 7.330071221881720772; Ej[8] = -0.276174333067751758;
            Kj[9] = 11.87151259742530180; Ej[9] = -0.393513114304375851;
            Kj[10] = 19.45851374822937738; Ej[10] = -0.575754406027879147;
            Kj[11] = 32.20638657246426863; Ej[11] = -0.860523235727239756;
            Kj[12] = 53.73749198700554656; Ej[12] = -1.308833205758540162;
            Kj[13] = 90.27388602940998849;
            for(j = 0; j <= 12; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.45;
            }
            *ellf += Kj[13] * powj;
            break;
        case 5:
            Kj[0] = 1.898924910271553526; Ej[0] = +1.325024497958230082;
            Kj[1] = 0.950521794618244435; Ej[1] = -0.521727647557566767;
            Kj[2] = 1.151077589959015808; Ej[2] = -0.194906430482126213;
            Kj[3] = 1.750239106986300540; Ej[3] = -0.171623726822011264;
            Kj[4] = 2.952676812636875180; Ej[4] = -0.202754652926419141;
            Kj[5] = 5.285800396121450889; Ej[5] = -0.278798953118534762;
            Kj[6] = 9.832485716659979747; Ej[6] = -0.420698457281005762;
            Kj[7] = 18.78714868327559562; Ej[7] = -0.675948400853106021;
            Kj[8] = 36.61468615273698145; Ej[8] = -1.136343121839229244;
            Kj[9] = 72.45292395127771801; Ej[9] = -1.976721143954398261;
            Kj[10] = 145.1079577347069102; Ej[10] = -3.531696773095722506;
            Kj[11] = 293.4786396308497026; Ej[11] = -6.446753640156048150;
            Kj[12] = 598.3851815055010179; Ej[12] = -11.97703130208884026;
            Kj[13] = 1228.420013075863451;
            Kj[14] = 2536.529755382764488;
            for(j = 0; j <= 12; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.55;
            }
            *ellf += Kj[13] * powj;
            powj *= m - 0.55;
            *ellf += Kj[14] * powj;
            break;
        case 6:
            Kj[0] = 2.007598398424376302; Ej[0] = +1.270707479650149744;
            Kj[1] = 1.248457231212347337; Ej[1] = -0.566839168287866583;
            Kj[2] = 1.926234657076479729; Ej[2] = -0.262160793432492598;
            Kj[3] = 3.751289640087587680; Ej[3] = -0.292244173533077419;
            Kj[4] = 8.119944554932045802; Ej[4] = -0.440397840850423189;
            Kj[5] = 18.66572130873555361; Ej[5] = -0.774947641381397458;
            Kj[6] = 44.60392484291437063; Ej[6] = -1.498870837987561088;
            Kj[7] = 109.5092054309498377; Ej[7] = -3.089708310445186667;
            Kj[8] = 274.2779548232413480; Ej[8] = -6.667595903381001064;
            Kj[9] = 697.5598008606326163; Ej[9] = -14.89436036517319078;
            Kj[10] = 1795.716014500247129; Ej[10] = -34.18120574251449024;
            Kj[11] = 4668.381716790389910; Ej[11] = -80.15895841905397306;
            Kj[12] = 12235.76246813664335; Ej[12] = -191.3489480762984920;
            Kj[13] = 32290.17809718320818; Ej[13] = -463.5938853480342030;
            Kj[14] = 85713.07608195964685; Ej[14] = -1137.380822169360061;
            Kj[15] = 228672.1890493117096;
            Kj[16] = 612757.2711915852774;
            for(j = 0; j <= 14; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.65;
            }
            *ellf += Kj[15] * powj;
            powj *= m - 0.65;
            *ellf += Kj[16] * powj;
//            printf("ellf elle %f %f\n",*ellf, *elle);
            break;
        case 7:
            Kj[0] = 2.156515647499643235; Ej[0] = +1.211056027568459525;
            Kj[1] = 1.791805641849463243; Ej[1] = -0.630306413287455807;
            Kj[2] = 3.826751287465713147; Ej[2] = -0.387166409520669145;
            Kj[3] = 10.38672468363797208; Ej[3] = -0.592278235311934603;
            Kj[4] = 31.40331405468070290; Ej[4] = -1.237555584513049844;
            Kj[5] = 100.9237039498695416; Ej[5] = -3.032056661745247199;
            Kj[6] = 337.3268282632272897; Ej[6] = -8.181688221573590762;
            Kj[7] = 1158.707930567827917; Ej[7] = -23.55507217389693250;
            Kj[8] = 4060.990742193632092; Ej[8] = -71.04099935893064956;
            Kj[9] = 14454.00184034344795; Ej[9] = -221.8796853192349888;
            Kj[10] = 52076.66107599404803; Ej[10] = -712.1364793277635425;
            Kj[11] = 189493.6591462156887; Ej[11] = -2336.125331440396407;
            Kj[12] = 695184.5762413896145; Ej[12] = -7801.945954775964673;
            Kj[13] = 2567994.048255284686; Ej[13] = -26448.19586059191933;
            Kj[14] = 9541921.966748386322; Ej[14] = -90799.48341621365251;
            Kj[15] = 35634927.44218076174; Ej[15] = -315126.0406449163424;
            Kj[16] = 133669298.4612040871; Ej[16] = -1104011.344311591159;
            Kj[17] = 503352186.6866284541;
            Kj[18] = 1901975729.538660119;
            Kj[19] = 7208915015.330103756;
            for(j = 0; j <= 16; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.75;
            }
            *ellf += Kj[17] * powj;
            powj *= m - 0.75;
            *ellf += Kj[18] * powj;
            powj *= m - 0.75;
            *ellf += Kj[19] * powj;
            break;
        case 8:
            Kj[0] = 2.318122621712510589; Ej[0] = +1.161307152196282836;
            Kj[1] = 2.616920150291232841; Ej[1] = -0.701100284555289548;
            Kj[2] = 7.897935075731355823; Ej[2] = -0.580551474465437362;
            Kj[3] = 30.50239715446672327; Ej[3] = -1.243693061077786614;
            Kj[4] = 131.4869365523528456; Ej[4] = -3.679383613496634879;
            Kj[5] = 602.9847637356491617; Ej[5] = -12.81590924337895775;
            Kj[6] = 2877.024617809972641; Ej[6] = -49.25672530759985272;
            Kj[7] = 14110.51991915180325; Ej[7] = -202.1818735434090269;
            Kj[8] = 70621.44088156540229; Ej[8] = -869.8602699308701437;
            Kj[9] = 358977.2665825309926; Ej[9] = -3877.005847313289571;
            Kj[10] = 1847238.263723971684; Ej[10] = -17761.70710170939814;
            Kj[11] = 9600515.416049214109; Ej[11] = -83182.69029154232061;
            Kj[12] = 50307677.08502366879; Ej[12] = -396650.4505013548170;
            Kj[13] = 265444188.6527127967; Ej[13] = -1920033.413682634405;
            Kj[14] = 1408862325.028702687;
            Kj[15] = 7515687935.373774627;
            for(j = 0; j <= 13; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.825;
            }
            *ellf += Kj[14] * powj;
            powj *= m - 0.825;
            *ellf += Kj[15] * powj;
            break;
        case 9:
            Kj[0] = 2.473596173751343912; Ej[0] = +1.124617325119752213;
            Kj[1] = 3.727624244118099310; Ej[1] = -0.770845056360909542;
            Kj[2] = 15.60739303554930496; Ej[2] = -0.844794053644911362;
            Kj[3] = 84.12850842805887747; Ej[3] = -2.490097309450394453;
            Kj[4] = 506.9818197040613935; Ej[4] = -10.23971741154384360;
            Kj[5] = 3252.277058145123644; Ej[5] = -49.74900546551479866;
            Kj[6] = 21713.24241957434256; Ej[6] = -267.0986675195705196;
            Kj[7] = 149037.0451890932766; Ej[7] = -1532.665883825229947;
            Kj[8] = 1043999.331089990839; Ej[8] = -9222.313478526091951;
            Kj[9] = 7427974.817042038995; Ej[9] = -57502.51612140314030;
            Kj[10] = 53503839.67558661151; Ej[10] = -368596.1167416106063;
            Kj[11] = 389249886.9948708474; Ej[11] = -2415611.088701091428;
            Kj[12] = 2855288351.100810619; Ej[12] = -16120097.81581656797;
            Kj[13] = 21090077038.76684053; Ej[13] = -109209938.5203089915;
            Kj[14] = 156699833947.7902014; Ej[14] = -749380758.1942496220;
            Kj[15] = 1170222242422.439893; Ej[15] = -5198725846.725541393;
            Kj[16] = 8777948323668.937971; Ej[16] = -36409256888.12139973;
            Kj[17] = 66101242752484.95041;
            Kj[18] = 499488053713388.7989;
            Kj[19] = 37859743397240299.20;
            for(j = 0; j <= 16; j++) {
                *ellf += Kj[j] * powj; 
                *elle += Ej[j] * powj; 
                powj *= m - 0.875;
            }
            *ellf += Kj[17] * powj;
            powj *= m - 0.875;
            *ellf += Kj[18] * powj;
            powj *= m - 0.875;
            *ellf += Kj[19] * powj;
            break;
        case 10:
            Kj[0] = 1.591003453790792180; Ej[0] = +1.550973351780472328;
            Kj[1] = 0.416000743991786912; Ej[1] = -0.400301020103198524;
            Kj[2] = 0.245791514264103415; Ej[2] = -0.078498619442941939;
            Kj[3] = 0.179481482914906162; Ej[3] = -0.034318853117591992;
            Kj[4] = 0.144556057087555150; Ej[4] = -0.019718043317365499;
            Kj[5] = 0.123200993312427711; Ej[5] = -0.013059507731993309;
            Kj[6] = 0.108938811574293531; Ej[6] = -0.009442372874146547;
            Kj[7] = 0.098853409871592910; Ej[7] = -0.007246728512402157;
            Kj[8] = 0.091439629201749751; Ej[8] = -0.005807424012956090;
            Kj[9] = 0.085842591595413900; Ej[9] = -0.004809187786009338;
            Kj[10] = 0.081541118718303215;
            for(j = 0; j <= 9; j++) {
                ellf1 += Kj[j] * powj; 
                powj *= 0.95 - m;
//                elle1 += Ej[j] * pow(0.95 - m, j); 
            }
            ellf1 += Kj[10] * powj;

            qj[1] = 1.0 / 16.0;
            qj[2] = 1.0 / 32.0;
            qj[3] = 21.0 / 1024.0;
            qj[4] = 31.0 / 2048.0;
            qj[5] = 6257.0 / 524288.0;
            qj[6] = 10293.0 / 1048576.0;
            qj[7] = 279025.0 / 33554432.0;
            qj[8] = 483127.0 / 67108864.0;
            qj[9] = 435506703.0 / 68719476736.0;
            qj[10] = 776957575.0 / 137438953472.0;
            qj[11] = 22417045555.0 / 4398046511104.0;
            qj[12] = 40784671953.0 / 8796093022208.0;
            qj[13] = 9569130097211.0 / 2251799813685248.0;
            qj[14] = 17652604545791.0 / 4503599627370496.0;

            powj = 1.0 - m;
            for(j = 1; j<= 14; j++) {
                q1 += qj[j] * powj;
                powj *= 1.0 - m;
            }

            Hj[0] = 0.040030102010319852;
            Hj[1] = 0.816301764094985436;
            Hj[2] = 0.324290133707045355;
            Hj[3] = 0.213800336032498154;
            Hj[4] = 0.164274100404920649;
            Hj[5] = 0.136260501044421020;
            Hj[6] = 0.118381184448440078;
            Hj[7] = 0.106100138383995067;
            Hj[8] = 0.097247053214705841;
            Hj[9] = 0.090651779381423238;
            Hj[10] = 0.085627517951558365;
           
            powj = 1.0;
            for(j = 0; j <= 10; j++) {
                Hm += Hj[j] * powj;
                powj *= 0.95 - m;
            }

            *ellf = -log(q1) * (ellf1/M_PI);
            *elle = 1.0 / ellf1 * (M_PI/2.0 + *ellf * Hm);
            break;
        default:
            printf("Error, bad input, quitting\n");
            break;
    }
}

/******************************************************************************
 * Calculate the magnetic and electric field at any point of the system.
 ******************************************************************************/
void getemf_wlcs(double x, double y, double z, double t, struct emfields *emf_tot)
{
    double B0; // background magnetic field
    struct emfields emf;
    struct bfields bmf;
    int i;
   
    emf.Bx = 0.0; emf.By = 0.0; emf.Bz = 0.0;
    emf.Ex = 0.0; emf.Ey = 0.0; emf.Ez = 0.0;
    bmf.Bx = 0.0; bmf.By = 0.0; bmf.Bz = 0.0;

    B0 = 0.0;
    for(i = 0; i < nwlcs; i++){
	    if(iBessel == 0){
		    getWireB(x, y, z, i, &bmf);
            emf_tot->Bx += bmf.Bx;
            emf_tot->By += bmf.By;
            emf_tot->Bz += bmf.Bz;
        }
	    else{
		    getWire_Bessel(x, y, z, i, t, &emf);
            emf_tot->Bx += emf.Bx; emf_tot->Ex += emf.Ex;
            emf_tot->By += emf.By; emf_tot->Ey += emf.Ey; 
            emf_tot->Bz += emf.Bz; emf_tot->Ez += emf.Ez;
        }

	    getLoopB(x, y, z, i, &bmf);

        emf_tot->Bx += bmf.Bx;
        emf_tot->By += bmf.By;
        emf_tot->Bz += bmf.Bz;
    }

    emf_tot->Bz += B0;
//    printf("%lf %lf %f %lf %lf %lf\n", emf_tot->Bx, emf_tot->By, emf_tot->Bz,
//            emf_tot->Ex, emf_tot->Ey, emf_tot->Ez);
}

/******************************************************************************
 * Transform coodinates between different reference frames.
 * Only two of three Euler angles are considered for the wire-loop system.                                         
 ******************************************************************************/
void transform(double *Ax, double *Ay, double *Az, 
        double cos_euler1, double sin_euler1, 
        double cos_euler2, double sin_euler2)
{
    double ax1, ay1, az1;
    
    ax1 =  cos_euler2 * (*Ax) + sin_euler2 * (*Ay);
    ay1 = -sin_euler2 * (*Ax) + cos_euler2 * (*Ay);
    az1 = *Az;

    *Ax = cos_euler1 * ax1 - sin_euler1 * az1;
    *Ay = ay1;
    *Az = sin_euler1 * ax1 + cos_euler1 * az1; 
}

/******************************************************************************
 * Transform fields between different reference frames. Only two of three Euler
 * angles are considered for the wire-loop system.                                         
 ******************************************************************************/
void reverse_tran(double *Ax, double *Ay, double *Az, 
        double cos_euler1, double sin_euler1, 
        double cos_euler2, double sin_euler2)
{
    double ax1, ay1, az1;
    
    ax1 = cos_euler1 * (*Ax) + sin_euler1 * (*Az);
    ay1 = *Ay;
    az1 = -sin_euler1 * (*Ax) + cos_euler1 * (*Az); 

    *Ax = cos_euler2 * ax1 - sin_euler2 * ay1;
    *Ay = sin_euler2 * ax1 + cos_euler2 * ay1;
    *Az = az1;
}
