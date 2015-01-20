#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "cuda.h"
#include "book.h"
#include "wb.h"
#include "Global.h"
#include <cuda_runtime.h>

#define DIM 256
//#define N (256*256*2)
#define N (256*64*3)
#define imin(a,b) (a<b?a:b)
const int threadsPerBlock = 256;
const int blocksPerGrid = imin( 256, (N/3+threadsPerBlock-1) / threadsPerBlock );

#define rnd( x ) (x * rand() / RAND_MAX)
#define INF 2e10f

const double ChargeMass = ChargeMass_proton;
const int nbins=110;            // number of energy bins 
const int bins_tot = nbins;     // two size of box
//const double center = 1.5;
//const double center = 0.75;
const double center = 1.25;
const double halfLen = 1.25;  // half of the side length of one box
const double Len = halfLen * 2.0;
const double Len2 = 1.5*Len;

__constant__ conf_sys config[Isys]; // constant memory

__global__ void kernel(int size, struct output_data *odat, struct output_data *InitData, 
                       double *DistSquare, int *PtclInBox, int *BinPtcl, int iBessel)
{
    unsigned int pn = threadIdx.x + blockIdx.x * blockDim.x;
    unsigned int i = threadIdx.x;
    //double dt = 1.0E-5;
    double dt = 1.0E-5; // for lower energy
    double epsilon = 1.0E-6;
    double tmp1 = 0.0;
    int tmp2 = 0;
    double vtot, ene;

    __shared__ double cache[threadsPerBlock];
    __shared__ int cache2[threadsPerBlock];
    __shared__ double x_ds[threadsPerBlock];
    __shared__ double y_ds[threadsPerBlock];
    __shared__ double z_ds[threadsPerBlock];
    __shared__ double vx_ds[threadsPerBlock];
    __shared__ double vy_ds[threadsPerBlock];
    __shared__ double vz_ds[threadsPerBlock];
    __shared__ double t_ds[threadsPerBlock];
    __shared__ double x_init[threadsPerBlock];
    __shared__ double y_init[threadsPerBlock];
    __shared__ double z_init[threadsPerBlock];
    __shared__ double t_init[threadsPerBlock];
    __shared__ int BinInOut[threadsPerBlock];
    __shared__ int In_Initial[threadsPerBlock];

    cache[i] = 0.0;
    cache2[i] = 0;
    In_Initial[i] = 0;
    __syncthreads();

    while(pn < size) {
        x_ds[i] = odat[pn].x;
        y_ds[i] = odat[pn].y;
        z_ds[i] = odat[pn].z;
        vx_ds[i] = odat[pn].vx;
        vy_ds[i] = odat[pn].vy;
        vz_ds[i] = odat[pn].vz;
        t_ds[i] = odat[pn].t;
        x_init[i] = InitData[pn].x;
        y_init[i] = InitData[pn].y;
        z_init[i] = InitData[pn].z;
        t_init[i] = InitData[pn].t;
        BinInOut[i] = BinPtcl[pn];

        __syncthreads();
        
        /*
        advance(&vx_ds[i], &vy_ds[i], &vz_ds[i],
                &x_ds[i], &y_ds[i], &z_ds[i], &t_ds[i], epsilon, dt, iBessel);
        In_Initial[i] = 1;

        __syncthreads();

        tmp1 = (pow(x_ds[i]-x_init[i], 2) + 
                pow(y_ds[i]-y_init[i], 2) + 
                pow(z_ds[i]-z_init[i], 2)) / (t_ds[i]-t_init[i]) / 2.0;
        tmp2 = 1;
        */
        // don't follow the particles getting out anymore
        if(fabs(x_ds[i]-center) < Len2 && fabs(y_ds[i]-center) < Len2 &&
           fabs(z_ds[i]-center) < Len2) {
            //advance(&vx_ds[i], &vy_ds[i], &vz_ds[i],
            //        &x_ds[i], &y_ds[i], &z_ds[i], &t_ds[i], epsilon, dt, iBessel);
            tracking_wirz(&vx_ds[i], &vy_ds[i], &vz_ds[i], iBessel,
                    &x_ds[i], &y_ds[i], &z_ds[i], &t_ds[i], dt);
            In_Initial[i] = 1;
        }
        __syncthreads();

        if(fabs(x_ds[i]-center) < Len2 && fabs(y_ds[i]-center) < Len2 &&
           fabs(z_ds[i]-center) < Len2) {
            tmp1 = ((x_ds[i]-x_init[i]) * (x_ds[i]-x_init[i]) +
                   (y_ds[i]-y_init[i]) * (y_ds[i]-y_init[i]) +
                   (z_ds[i]-z_init[i]) * (z_ds[i]-z_init[i])) / ((t_ds[i]-t_init[i]) * 2);
            tmp2 = 1;
        }

        cache[i] += tmp1;
        cache2[i] += tmp2;

        __syncthreads();
        
        if(In_Initial[i] == 1) {
            vtot = vx_ds[i]*vx_ds[i] + vy_ds[i]*vy_ds[i] + vz_ds[i]*vz_ds[i];
            ene = (1.0/sqrt(1.0-vtot) - 1.0) * pene;
            BinInOut[i] = floor((log10(ene) + 7.0) / 0.1);
        }

        __syncthreads();

        odat[pn].t = t_ds[i];
        odat[pn].x = x_ds[i];
        odat[pn].y = y_ds[i];
        odat[pn].z = z_ds[i];
        odat[pn].vx = vx_ds[i];
        odat[pn].vy = vy_ds[i];
        odat[pn].vz = vz_ds[i];
        BinPtcl[pn] = BinInOut[i];

        pn += blockDim.x * gridDim.x;
    }

    int stride = blockDim.x/2;

    while (stride != 0) {
        __syncthreads();
        if(i < stride) {
            cache[i] += cache[i + stride];
            cache2[i] += cache2[i + stride];
        }
        stride /= 2;
    }
    __syncthreads();

    if (i == 0) {
        DistSquare[blockIdx.x] = cache[0];
        PtclInBox[blockIdx.x] = cache2[0];
    }
};

/******************************************************************************
 * Kernel to sort the simulation result to avoid control divergence.
 * Reference: https://gist.github.com/1392067
 ******************************************************************************/
__global__ void bitonic_sort_step(struct output_data *odat, struct output_data *InitData, 
                                  int size, int j, int k, int *BinPtcl)
{
    unsigned int i, ixj;
    i = threadIdx.x + blockIdx.x * blockDim.x;
    ixj = i^j;
    
    while (i < size) {
        if((ixj)>i) {
            if((i&k)==0) {
                /* Sort ascending */
                if((fabs(odat[i].x-center) > Len2 || 
                    fabs(odat[i].y-center) > Len2 || 
                    fabs(odat[i].z-center) > Len2) &&
                   (fabs(odat[ixj].x-center) < Len2 && 
                    fabs(odat[ixj].y-center) < Len2 && 
                    fabs(odat[ixj].z-center) < Len2)) {
                    /* exchange(i,ixj) */
                    swap_double(&odat[i].t, &odat[ixj].t);
                    swap_double(&odat[i].x, &odat[ixj].x);
                    swap_double(&odat[i].y, &odat[ixj].y);
                    swap_double(&odat[i].z, &odat[ixj].z);
                    swap_double(&odat[i].vx, &odat[ixj].vx);
                    swap_double(&odat[i].vy, &odat[ixj].vy);
                    swap_double(&odat[i].vz, &odat[ixj].vz);
                    swap_double(&InitData[i].t, &InitData[ixj].t);
                    swap_double(&InitData[i].x, &InitData[ixj].x);
                    swap_double(&InitData[i].y, &InitData[ixj].y);
                    swap_double(&InitData[i].z, &InitData[ixj].z);
                    swap_double(&InitData[i].vx, &InitData[ixj].vx);
                    swap_double(&InitData[i].vy, &InitData[ixj].vy);
                    swap_double(&InitData[i].vz, &InitData[ixj].vz);
                    swap_int(&BinPtcl[i], &BinPtcl[ixj]);
                }
            }

            if((i&k)!=0) {
                /* Sort descending */
                if((fabs(odat[ixj].x-center) > Len2 || 
                    fabs(odat[ixj].y-center) > Len2 || 
                    fabs(odat[ixj].z-center) > Len2) &&
                   (fabs(odat[i].x-center) < Len2 && 
                    fabs(odat[i].y-center) < Len2 && 
                    fabs(odat[i].z-center) < Len2)) {
                    /* exchange(i,ixj) */
                    swap_double(&odat[i].t, &odat[ixj].t);
                    swap_double(&odat[i].x, &odat[ixj].x);
                    swap_double(&odat[i].y, &odat[ixj].y);
                    swap_double(&odat[i].z, &odat[ixj].z);
                    swap_double(&odat[i].vx, &odat[ixj].vx);
                    swap_double(&odat[i].vy, &odat[ixj].vy);
                    swap_double(&odat[i].vz, &odat[ixj].vz);
                    swap_double(&InitData[i].t, &InitData[ixj].t);
                    swap_double(&InitData[i].x, &InitData[ixj].x);
                    swap_double(&InitData[i].y, &InitData[ixj].y);
                    swap_double(&InitData[i].z, &InitData[ixj].z);
                    swap_double(&InitData[i].vx, &InitData[ixj].vx);
                    swap_double(&InitData[i].vy, &InitData[ixj].vy);
                    swap_double(&InitData[i].vz, &InitData[ixj].vz);
                    swap_int(&BinPtcl[i], &BinPtcl[ixj]);
                }
            }
        }
        i += blockDim.x * gridDim.x; 
    }
};

int main(int argc, char ** argv)
{
    wbArg_t args;
    int deviceCount;
    HANDLE_ERROR(cudaGetDeviceCount(&deviceCount));
    if (deviceCount < 2) {
        printf( "We need at least two compute 1.0 or greater "
                "devices, but only found %d\n", deviceCount );
        return 0;
    }

    args = wbArg_read(argc, argv);
    wbTime_start(Generic, "Importing data and creating memory on host");
    
    int iBessel;
    // allocate memory for the initial data
	conf_sys *config_tmp = (conf_sys*)malloc(sizeof(conf_sys)*Isys);
	ptl_init *init = (ptl_init*)malloc(sizeof(ptl_init));
	read_init(config_tmp, init, &iBessel);

    srand48(time(NULL));
    double phi, theta;
    double randn; // normally distributed random number
    double sigma;
    
    DataStruct data[deviceCount];
    for(int i = 0; i < deviceCount; i++) {
        data[i].size = N/deviceCount;
    }
    // Create streams for issuing GPU command asynchronously and 
    // allocate memory (GPU and System page-locked)
    double vtot, ene;
    for(int i = 0; i < deviceCount; i++) {
        HANDLE_ERROR(cudaSetDevice(i));
        HANDLE_ERROR(cudaStreamCreate(&data[i].stream));
        // Allocate memory
        HANDLE_ERROR(cudaMalloc((void**)&data[i].output_d, sizeof(output_data)*data[i].size));
        HANDLE_ERROR(cudaMalloc((void**)&data[i].InitData, sizeof(output_data)*data[i].size));
        HANDLE_ERROR(cudaMalloc((void**)&data[i].DistSquare_d, sizeof(double)* blocksPerGrid));
        HANDLE_ERROR(cudaMalloc((void**)&data[i].PtclInBox_d, sizeof(int)* blocksPerGrid));
        HANDLE_ERROR(cudaMalloc((void**)&data[i].BinPtcl_d, sizeof(int)*data[i].size));
        HANDLE_ERROR(cudaMallocHost((void**)&data[i].BinPtcl_h, sizeof(int) * data[i].size));
        HANDLE_ERROR(cudaMallocHost((void**)&data[i].BinInit, sizeof(int) * data[i].size));
        HANDLE_ERROR(cudaMallocHost((void**)&data[i].output_h, sizeof(output_data) * data[i].size));
        HANDLE_ERROR(cudaMallocHost((void**)&data[i].DistSquare_h, sizeof(double) * blocksPerGrid));
        HANDLE_ERROR(cudaMallocHost((void**)&data[i].PtclInBox_h, sizeof(int) * blocksPerGrid));
        //double v0 = 3.2646E-4;
        double v0 = 7.3E-4;
        for(int j = 0; j < data[i].size; j++) {
            data[i].output_h[j].t = 0.0;
            data[i].output_h[j].x = drand48() * Len;
            data[i].output_h[j].y = drand48() * Len;
            data[i].output_h[j].z = drand48() * Len;
            /*
            theta = drand48() * M_PI;
            phi = drand48() * 2.0 * M_PI;
            
            data[i].output_h[j].vx = v0 * sin(theta) * cos(phi);
            data[i].output_h[j].vy = v0 * sin(theta) * sin(phi);
            data[i].output_h[j].vz = v0 * cos(theta);
            */
            // Maxwellian distribution of velocity for each components
            sigma = 9.09E4; // sigma = sqrt(kT/m)
            randn = erfinv(2.0*drand48() - 1.0) * sqrt(2.0) * sigma;
            data[i].output_h[j].vx = randn * 100 / c0; // CGS
            randn = erfinv(2.0*drand48() - 1.0) * sqrt(2.0) * sigma;
            data[i].output_h[j].vy = randn * 100 / c0; // CGS
            randn = erfinv(2.0*drand48() - 1.0) * sqrt(2.0) * sigma;
            data[i].output_h[j].vz = randn * 100 / c0; // CGS

            vtot = data[i].output_h[j].vx * data[i].output_h[j].vx + 
                   data[i].output_h[j].vy * data[i].output_h[j].vy + 
                   data[i].output_h[j].vz * data[i].output_h[j].vz;
            ene = (1.0/sqrt(1.0-vtot) - 1.0) * pene;
            data[i].BinInit[j] = floor((log10(ene)+7.0)/0.1);
        }
    }
    
    //Copy data to GPU, launch the kernel and copy data back. All asynchronously
    wbTime_start(Compute, "Performing CUDA computation");

    for(int i = 0; i < deviceCount; i++) {
        HANDLE_ERROR(cudaSetDevice(i));
        //Copy input data from CPU
        HANDLE_ERROR(cudaMemcpyAsync(data[i].output_d, data[i].output_h, data[i].size*sizeof(output_data), 
                                    cudaMemcpyHostToDevice, data[i].stream));
        HANDLE_ERROR(cudaMemcpyAsync(data[i].InitData, data[i].output_h, data[i].size*sizeof(output_data), 
                                    cudaMemcpyHostToDevice, data[i].stream));
        HANDLE_ERROR(cudaMemcpyAsync(data[i].BinPtcl_d, data[i].BinInit, data[i].size*sizeof(int), 
                                     cudaMemcpyHostToDevice, data[i].stream));
        HANDLE_ERROR(cudaMemcpyToSymbol(config, config_tmp, sizeof(conf_sys)*Isys));
    }

    FILE *fp2, *fpIn, *fpOut;
    char buf1[256], buf2[256], buf3[256];
    snprintf(buf1, sizeof buf1, "%s", "drr.txt");
    snprintf(buf2, sizeof buf2, "%s", "SpectIn.txt");
    snprintf(buf3, sizeof buf3, "%s", "SpectOut.txt");
    fp2 = fopen(buf1, "w");
    fpIn = fopen(buf2, "w");
    fpOut = fopen(buf3, "w");
    fprintf(fp2, "%d %10.8e\n", N, 0.0);
   
    int *spectIn, *spectOut;
    spectIn = (int*)malloc(sizeof(int)*bins_tot);
    spectOut = (int*)malloc(sizeof(int)*bins_tot);
    for(int i = 0; i < bins_tot; i++) {
        spectIn[i] = 0;
        spectOut[i] = 0;
    }
    for(int i = 0; i < deviceCount; i++) {
        for(int j = 0; j < data[i].size; j++) {
            spectIn[data[i].BinInit[j]]++;
        }
    }
    for(int i = 0; i < bins_tot; i++) {
        fprintf(fpIn, "%d ", spectIn[i]);
        fprintf(fpOut, "%d ", spectOut[i]);
    }
    fprintf(fpIn, "\n");
    fprintf(fpOut, "\n");

    int step = 0;
    int bj, bk;
    while(step < 2E6) {
        //printf("%d\n", step);
        //Perform GPU computation
        for(int i = 0; i < deviceCount; i++) {
            HANDLE_ERROR(cudaSetDevice(i));
            kernel<<<blocksPerGrid,threadsPerBlock,0,data[i].stream>>>(data[i].size, data[i].output_d, 
                                                     data[i].InitData, data[i].DistSquare_d, data[i].PtclInBox_d,
                                                     data[i].BinPtcl_d, iBessel);
            cudaThreadSynchronize();
        }
        step += 1;
        if(step % 10 == 0) {
            for(int i = 0; i < deviceCount; i++) {
                HANDLE_ERROR(cudaSetDevice(i));
                for(bk = 2; bk <= data[i].size; bk<<=1) {
                    for(bj=bk>>1; bj > 0; bj=bj>>1) {
                        bitonic_sort_step<<<blocksPerGrid, threadsPerBlock,0,data[i].stream>>>
                                        (data[i].output_d, data[i].InitData, data[i].size, bj, bk, data[i].BinPtcl_d);
                    }
                }
                HANDLE_ERROR(cudaMemcpyAsync(data[i].DistSquare_h, data[i].DistSquare_d, blocksPerGrid*sizeof(double), 
                                            cudaMemcpyDeviceToHost, data[i].stream));
                HANDLE_ERROR(cudaMemcpyAsync(data[i].PtclInBox_h, data[i].PtclInBox_d, blocksPerGrid*sizeof(int), 
                                            cudaMemcpyDeviceToHost, data[i].stream));
            }
            for(int i = 0; i < deviceCount; i++) {
                HANDLE_ERROR(cudaSetDevice(i));
                //Wait for all operations to finish
                cudaStreamSynchronize(data[i].stream);
                double drr = 0.0;
                int PtclN = 0;
                for(int j = 0; j < blocksPerGrid; j++) {
                    drr += data[i].DistSquare_h[j];
                    PtclN += data[i].PtclInBox_h[j];
                    //printf("%10.8e\n", data[i].DistSquare_h[j]);
                }
                data[i].drr = drr;
                data[i].PtclN = PtclN;
            }
            double sumDrr = 0.0;
            int PtclTot = 0;
            for(int i = 0; i < deviceCount; i++) {
                sumDrr += data[i].drr;
                PtclTot += data[i].PtclN;
            }
            fprintf(fp2, "%d %10.8e\n", PtclTot, sumDrr);
        }

        if(step % 100 == 0) {
            for(int i = 0; i < deviceCount; i++) {
                HANDLE_ERROR(cudaSetDevice(i));
                cudaThreadSynchronize();
                HANDLE_ERROR(cudaMemcpyAsync(data[i].BinPtcl_h, data[i].BinPtcl_d, data[i].size*sizeof(int), 
                                            cudaMemcpyDeviceToHost, data[i].stream));
            }
            // spectrum
            for (int i = 0; i < bins_tot; i++) {
                spectIn[i] = 0;
                spectOut[i] = 0;
            }
            for(int i = 0; i < deviceCount; i++) {
                for(int j = 0; j < data[i].PtclN; j++) {
                    spectIn[data[i].BinPtcl_h[j]]++;
                }
                for(int j = data[i].PtclN; j < data[i].size; j++) {
                    spectOut[data[i].BinPtcl_h[j]]++;
                }
            }
            for (int i = 0; i < bins_tot; i++) {
                fprintf(fpIn, "%d ", spectIn[i]);
                fprintf(fpOut, "%d ", spectOut[i]);
            }
            fprintf(fpIn, "\n");
            fprintf(fpOut, "\n");
        }
    }
    fclose(fp2);
    fclose(fpIn);
    fclose(fpOut);
    wbTime_stop(Compute, "Performing CUDA computation");

    free(config_tmp);
    free(spectIn);
    free(spectOut);
    for(int i = 0; i < deviceCount; i++) {
        HANDLE_ERROR(cudaSetDevice(i));
        HANDLE_ERROR(cudaFreeHost(data[i].output_h));
        HANDLE_ERROR(cudaFreeHost(data[i].DistSquare_h));
        HANDLE_ERROR(cudaFreeHost(data[i].PtclInBox_h));
        HANDLE_ERROR(cudaFreeHost(data[i].BinPtcl_h));
        HANDLE_ERROR(cudaFreeHost(data[i].BinInit));
        HANDLE_ERROR(cudaFree(data[i].output_d));
        HANDLE_ERROR(cudaFree(data[i].InitData));
        HANDLE_ERROR(cudaFree(data[i].DistSquare_d));
        HANDLE_ERROR(cudaFree(data[i].PtclInBox_d));
        HANDLE_ERROR(cudaFree(data[i].BinPtcl_d));
        HANDLE_ERROR(cudaStreamDestroy(data[i].stream));
    }
	return 0;
}

/******************************************************************************
 * Read the configuration data of a wire-loop current system (WLCS).
 * Sin and Cos functions are calculated here since they  will be re-used. 
 ******************************************************************************/
void read_wlcs(struct conf_sys *config, struct ptl_init *init, int *iBessel)
{
	FILE *fp;
    double theta, phi, alpha, beta;
    double omega;
    int i;
	//fp = fopen("init6.dat","r");
	fp = fopen("init8_3l.dat","r");
	//fp = fopen("init8_3l_001.dat","r");
	//fp = fopen("init8_3l_00001Hz.dat","r");
	//fp = fopen("init3_2.dat","r");
	for(i = 0; i<Isys; i++) {
		fscanf(fp, "%lf", &(config[i].cur_wr));
		fscanf(fp, "%lf", &(config[i].x_wr));
		fscanf(fp, "%lf", &(config[i].y_wr));
		fscanf(fp, "%lf", &(config[i].z_wr));
		fscanf(fp, "%lf", &theta);
		fscanf(fp, "%lf", &phi);
		fscanf(fp, "%lf", &(config[i].t0));
		fscanf(fp, "%lf", &(config[i].cur_lp));
		fscanf(fp, "%lf", &(config[i].r_lp));
		fscanf(fp, "%lf", &(config[i].x_lp));
		fscanf(fp, "%lf", &(config[i].y_lp));
		fscanf(fp, "%lf", &(config[i].z_lp));
		fscanf(fp, "%lf", &alpha);
		fscanf(fp, "%lf", &beta);

        // avoiding re-calculation for magntic field
		config[i].costheta_wr = cos(theta);
		config[i].sintheta_wr = sin(theta);
		config[i].cosphi_wr = cos(phi);
		config[i].sinphi_wr = sin(phi);
		config[i].cosalpha_lp = cos(alpha);
		config[i].sinalpha_lp = sin(alpha);
		config[i].cosbeta_lp = cos(beta);
		config[i].sinbeta_lp = sin(beta);
	}

	fscanf(fp, "%lf", &(init->x0));
	fscanf(fp, "%lf", &(init->y0));
	fscanf(fp, "%lf", &(init->z0));
	fscanf(fp, "%d", &(init->Nturn));
	fscanf(fp, "%d", &(init->Nstep));
	fscanf(fp, "%d", iBessel);
	fscanf(fp, "%lf", &(init->vsp));
	fscanf(fp, "%lf", &(init->t));
	fscanf(fp, "%lf", &(omega));

    for(i = 0; i < Isys; i++) {
        config[i].omega = omega;
    }

    /*
    for(i = 0; i < Isys; i++){
		printf("%lf ", config[i].cur_wr);
		printf("%lf ", config[i].x_wr);
		printf("%lf ", config[i].y_wr);
		printf("%lf ", config[i].z_wr);
		printf("%lf ", config[i].sintheta_wr);
		printf("%lf ", config[i].costheta_wr);
		printf("%lf ", config[i].sinphi_wr);
		printf("%lf ", config[i].cosphi_wr);
		printf("%lf ", config[i].t0);
		printf("%lf ", config[i].cur_lp);
		printf("%lf ", config[i].r_lp);
		printf("%lf ", config[i].x_lp);
		printf("%lf ", config[i].y_lp);
		printf("%lf ", config[i].z_lp);
		printf("%lf ", config[i].sinalpha_lp);
		printf("%lf ", config[i].cosalpha_lp);
		printf("%lf ", config[i].sinbeta_lp);
		printf("%lf ", config[i].cosbeta_lp);
        printf("%lf ", config[i].omega);
        printf("\n");
    }
    */
	fclose(fp);
}

/******************************************************************************
 * Calculate the magnetic field of one straight wire in time-dependent condition.
 ******************************************************************************/
void getWire_Bessel(double x, double y, double z, int iconf, double t,
        double *Bx, double *By, double *Bz, double *Ex, double *Ey, double *Ez)
{
	double xp, yp, zp;
	double x0, rho, const0;
	double bj0, bj1, by0, by1; // Bessel functions
	double wireB, wireE;
	double phi1, curI, omega;
    double t1;

    *Bx = 0.0; *By = 0.0; *Bz = 0.0;
    *Ex = 0.0; *Ey = 0.0; *Ez = 0.0;

	xp = x - config[iconf].x_wr;
	yp = y - config[iconf].y_wr;
	zp = z - config[iconf].z_wr;
    t1 = config[iconf].t0 + t;
    curI = config[iconf].cur_wr;
    omega = config[iconf].omega;

    transform(&xp, &yp, &zp, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr,
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
	rho = sqrt(xp*xp + yp*yp);   // distance
    x0 = 2.32E-2 * rho * omega; // page 9
/*
	bj0 = bessj0(x0);
	by0 = bessy0(x0);
	bj1 = bessj1(x0);
	by1 = bessy1(x0);
*/
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
		
	*Bx = wireB * cos(phi1);
	*By = wireB * sin(phi1);
	*Bz = 0.0;
		
    reverse_tran(Bx, By, Bz, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr, 
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);

//---------------------------------- Ex, Ey, Ez ----------------------------------
	wireE = -const0 * x0 * (by0*sin(omega * t1) + bj0*cos(omega * t1));

	*Ex = 0.0;
	*Ey = 0.0;
	*Ez = wireE;

    reverse_tran(Ex, Ey, Ez, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr,
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
}

/*-------------------------------------------------------
 * Calculate the magnetic field of one straight wire in
 * time-independent condition.
 *-------------------------------------------------------*/
void getWireB(double x, double y, double z, int iconf,
        double *Bx, double *By, double *Bz)
{
	double xp, yp, zp;
	double rho, cos_phi0, phi0;
	double Bphi, curI;

    *Bx = 0.0; *By = 0.0; *Bz = 0.0; // initilization
	
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

	*Bx = Bphi * cos(phi0 + M_PI/2.0);
	*By = Bphi * sin(phi0 + M_PI/2.0);
	*Bz = 0.0;

    reverse_tran(Bx, By, Bz, 
            config[iconf].costheta_wr, config[iconf].sintheta_wr, 
            config[iconf].cosphi_wr, config[iconf].sinphi_wr);
}

/******************************************************************************
 * Calculate the magnetic field of one loop in time-independent condition.
 ******************************************************************************/
void getLoopB(double x, double y, double z, int iconf,
        double *Bx, double *By, double *Bz)
{
	double xp, yp, zp, curI, al;
	double k_tmp, E_tmp, F_tmp;
	double x_tmp, rho_tmp, k2_tmp, BpRho4, BpRho3;
    double tmp1, tmp2, tmp3, tmp4;
    double rhop1, rhop3;
    double cosalpha, sinalpha, cosbeta, sinbeta;
	
    *Bx = 0.0; *By = 0.0; *Bz = 0.0; // initilization

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
	*Bx = BpRho3 * xp;
	*By = BpRho3 * yp;
	*Bz = tmp1 * (F_tmp - (tmp2 - tmp3) * tmp4);
        
    reverse_tran(Bx, By, Bz, cosbeta, sinbeta, cosalpha, sinalpha);
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
 * Some functions used in code, including cross product of two vectors and the
 * modulus of a vector.
 ******************************************************************************/
void cross_product(double beta_x, double beta_y, double beta_z,
		double Bx, double By, double Bz,
		double* Cx, double* Cy, double* Cz)
{
	*Cx = beta_y*Bz - beta_z*By;
	*Cy = beta_z*Bx - beta_x*Bz;
	*Cz = beta_x*By - beta_y*Bx;
}

double modulus(double x, double y, double z)
{
    return sqrt(x*x+y*y+z*z);
};

/******************************************************************************
 * Calculate the magnetic and electric field at any point of the system.
 ******************************************************************************/
void getEM(double x, double y, double z, double t, int iBessel,
        double *B_totx, double *B_toty, double *B_totz,
        double *E_totx, double *E_toty, double *E_totz)
{
    double Bx, By, Bz, Ex, Ey, Ez;
    double B0; // background magnetic field
    int i;
    
    B0 = 0.0;
    *B_totx = 0.0; *B_toty = 0.0; *B_totz = 0.0;
    *E_totx = 0.0; *E_toty = 0.0; *E_totz = 0.0;
    for(i = 0; i < Isys; i++){
	    if(iBessel == 0){
		    getWireB(x, y, z, i, &Bx, &By, &Bz);
            *B_totx += Bx;
            *B_toty += By;
            *B_totz += Bz;
        }
	    else{
		    getWire_Bessel(x, y, z, i, t, &Bx, &By, &Bz, &Ex, &Ey, &Ez);
            *B_totx += Bx; *E_totx += Ex;
            *B_toty += By; *E_toty += Ey; 
            *B_totz += Bz; *E_totz += Ez;
        }

	    getLoopB(x, y, z, i, &Bx, &By, &Bz);

        *B_totx += Bx;
        *B_toty += By;
        *B_totz += Bz;
    }
/*
    // force-free field as background
    const double A = 1.0, B = sqrt(2.0/3.0), C = sqrt(1.0/3.0);
    const double lamda = 1.0;
    *B_totx += A*sin(lamda*z) + C*cos(lamda*y);
    *B_toty += B*sin(lamda*x) + A*cos(lamda*z);
    *B_totz += C*sin(lamda*y) + B*cos(lamda*x);
*/

    *B_totz += B0;
}


/******************************************************************************
 * ============== Particle Tracking based on Wirz's method ==============
 * Reference:
 * Comparison of Charged Particle Tracking Methods for Non-Uniform Magnetic Fields,
 * Hann-Shin Mao and Richard E. Wirz, June 2011
 *
 * Wirz, R., Discharge plasma processes of ring-cusp ion thrusters, 
 * Ph.D. Thesis, California Institute of Technology, 2005.
 *
 * Author: Xiaocan Li
 * Date: Oct-26-2012
 ******************************************************************************/

void tracking_wirz(double *vx, double *vy, double *vz, int iBessel, 
                              double *x, double *y, double *z, double *t, double dt)
{
//    double epsilon = 1.0E-7;

    double h, B_tot;
    double B_totx, B_toty, B_totz;
    double E_totx, E_toty, E_totz;
    double vx_minus, vx_plus, vxold;
    double vy_minus, vy_plus, vyold;
    double vz_minus, vz_plus, vzold;
    double deltax_m, deltay_m, deltaz_m;
    double deltax_new, deltay_new, deltaz_new;
    double hparax, hparay, hparaz; // prallel unit vector to B
    double hperpx, hperpy, hperpz; // perpendicular unit vector
    double hrx, hry, hrz; // radial unit vector
    double dx_para, dy_para, dz_para, dx_perp, dy_perp, dz_perp, dx_r, dy_r, dz_r;
    double dtheta, gyroR, omegac; // gyro frequency and Larmor radius
    double xpre, ypre, zpre; // predicted midpoint
    double vperp; // perpendicular velocity
    double vdothp; // dot product of v and hp
    double x0, y0, z0, t0;
    double deltax, deltay, deltaz;
    double gama;

    x0 = *x; y0 = *y; z0 = *z; t0 = *t;
    vxold = *vx; vyold = *vy; vzold = *vz;
	getdelta(&deltax, &deltay, &deltaz, vxold, vyold, vzold);	
	gama = getgamma(deltax, deltay, deltaz);
//
    getEM(x0, y0, z0, t0, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);
//    B_totx = 0.0001; B_toty = 0.0003; B_totz = 0.00000001;
    B_tot = modulus(B_totx, B_toty, B_totz);
//    dt = epsilon / B_tot; // change dt for each position dt = epsilon * m / (q*B)
    h = dt / 2.0;
//    printf("%19.18e\n", dt);
//    printf("%10.9e\n", E_totz);

// calculate the first intermidate velocity v_minus by applying the first electric half impulse
    deltax_m = deltax + ChargeMass * h * E_totx;
    deltay_m = deltay + ChargeMass * h * E_toty;
    deltaz_m = deltaz + ChargeMass * h * E_totz;

//    printf("vold %20.19e\n ", vxold * vxold + vyold * vyold + vzold *vzold);

//    printf("deltax1 deltay1 deltaz1 %20.19e %20.19e %20.19e\n ", deltax, deltay, deltaz);
//    printf("deltax_m deltay_m deltaz_m %20.19e %20.19e %20.19e\n ", deltax_m, deltay_m, deltaz_m);

	gama = getgamma(deltax, deltay, deltaz);
    vxold = deltax / gama;
    vyold = deltay / gama;
    vzold = deltaz / gama;
//    printf("vold %20.19e\n ", vxold * vxold + vyold * vyold + vzold *vzold);
	gama = getgamma(deltax_m, deltay_m, deltaz_m);
    vx_minus = deltax_m / gama;
    vy_minus = deltay_m / gama;
    vz_minus = deltaz_m / gama;

// vminus is used to calculate a predicted midpoint
    hparax = B_totx / B_tot;
    hparay = B_toty / B_tot;
    hparaz = B_totz / B_tot;
    vdothp = vx_minus * hparax + vy_minus * hparay + vz_minus * hparaz; 
    vperp = modulus(vx_minus - vdothp * hparax, vy_minus - vdothp * hparay, 
                    vz_minus - vdothp * hparaz);
    hperpx = (vx_minus - vdothp * hparax) / vperp;  
    hperpy = (vy_minus - vdothp * hparay) / vperp;  
    hperpz = (vz_minus - vdothp * hparaz) / vperp; 
// for postive charge only, if it is negative charge, there should be minus there
    cross_product(hparax, hparay, hparaz, hperpx, hperpy, hperpz, &hrx, &hry, &hrz);

    gyroR = 1.0 / ChargeMass * vperp * gama * c0 / L0 / B_tot; // Larmor radius
    omegac = ChargeMass * B_tot / gama; // gyro frequency
    dtheta = omegac * dt;
    dx_para = h * vdothp * c0 * hparax / L0;
    dy_para = h * vdothp * c0 * hparay / L0;
    dz_para = h * vdothp * c0 * hparaz / L0;

    dx_perp = gyroR * sin(dtheta / 2.0) * hperpx;
    dy_perp = gyroR * sin(dtheta / 2.0) * hperpy;
    dz_perp = gyroR * sin(dtheta / 2.0) * hperpz;

    dx_r = -gyroR * (1.0 - cos(dtheta / 2.0)) * hrx;
    dy_r = -gyroR * (1.0 - cos(dtheta / 2.0)) * hry;
    dz_r = -gyroR * (1.0 - cos(dtheta / 2.0)) * hrz;

// predicted midpoint
    xpre = x0 + dx_para + dx_perp + dx_r;
    ypre = y0 + dy_para + dy_perp + dy_r;
    zpre = z0 + dz_para + dz_perp + dz_r;

// updata the coodinate using the predicted midpoint
//    B_totx = 0.0001; B_toty = 0.0003; B_totz = 0.00000001;
    getEM(xpre, ypre, zpre, t0 + h, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);
    B_tot = modulus(B_totx, B_toty, B_totz);

// vminus is used to calculate a predicted midpoint
    hparax = B_totx / B_tot;
    hparay = B_toty / B_tot;
    hparaz = B_totz / B_tot;
    
    vdothp = vx_minus * hparax + vy_minus * hparay + vz_minus * hparaz; 
    vperp = modulus(vx_minus - vdothp * hparax, vy_minus - vdothp * hparay, 
                    vz_minus - vdothp * hparaz);
    hperpx = (vx_minus - vdothp * hparax) / vperp;  
    hperpy = (vy_minus - vdothp * hparay) / vperp;  
    hperpz = (vz_minus - vdothp * hparaz) / vperp; 

// for postive charge only, if it is negative charge, there should be minus there
    cross_product(hparax, hparay, hparaz, hperpx, hperpy, hperpz, &hrx, &hry, &hrz);
//    printf("hr %10.9e\n ", hperpx * hrx +  hperpy * hry + hperpz * hrz);
    omegac = ChargeMass * B_tot / gama; // gyro frequency
    dtheta = omegac * dt;
    
    vx_plus = hparax * vdothp + vperp * (cos(dtheta) * hperpx - sin(dtheta) * hrx);
    vy_plus = hparay * vdothp + vperp * (cos(dtheta) * hperpy - sin(dtheta) * hry);
    vz_plus = hparaz * vdothp + vperp * (cos(dtheta) * hperpz - sin(dtheta) * hrz);
//
// update new velocity
// 
	getdelta(&deltax, &deltay, &deltaz, vx_plus, vy_plus, vz_plus);	
    deltax_new = deltax + ChargeMass * h * E_totx;
    deltay_new = deltay + ChargeMass * h * E_toty;
    deltaz_new = deltaz + ChargeMass * h * E_totz;

	gama = getgamma(deltax_new, deltay_new, deltaz_new);
    *vx = deltax_new / gama;
    *vy = deltay_new / gama;
    *vz = deltaz_new / gama;

//    printf("%20.19e %20.19e %20.19e\n ", vxold * vxold + vyold * vyold + vzold *vzold,
//            vx_minus * vx_minus + vy_minus * vy_minus + vz_minus * vz_minus,
//            *vx * *vx + *vy * *vy + *vz * *vz);
//
// update new position
//
    vdothp = vx_plus * hparax + vy_plus * hparay + vz_plus * hparaz; 
    vperp = modulus(vx_plus - vdothp * hparax, vy_plus - vdothp * hparay, 
                    vz_plus - vdothp * hparaz);
    gyroR = 1.0 / ChargeMass * vperp * gama * c0 / L0 / B_tot; // Larmor radius
    dx_para = dt * vdothp * c0 * hparax / L0;
    dy_para = dt * vdothp * c0 * hparay / L0;
    dz_para = dt * vdothp * c0 * hparaz / L0;

    dx_perp = gyroR * sin(dtheta) * hperpx;
    dy_perp = gyroR * sin(dtheta) * hperpy;
    dz_perp = gyroR * sin(dtheta) * hperpz;

    dx_r = -gyroR * (1.0 - cos(dtheta)) * hrx;
    dy_r = -gyroR * (1.0 - cos(dtheta)) * hry;
    dz_r = -gyroR * (1.0 - cos(dtheta)) * hrz;

    *x += dx_para + dx_perp + dx_r;
    *y += dy_para + dy_perp + dy_r;
    *z += dz_para + dz_perp + dz_r;
    *t += dt;
}

/******************************************************************************
 * 4th order Runge-Kutta method to track particles
 ******************************************************************************/
void advance(double *vx, double *vy, double *vz, double *x, double *y, double *z, 
		                double *t, double epsilon, double dt, int iBessel)
{
	double h;
	double B_totx, B_toty, B_totz;
	double E_totx, E_toty, E_totz;
	double vx0, vx1, vx2, vx3;
	double vy0, vy1, vy2, vy3;
	double vz0, vz1, vz2, vz3;
	double x0, x1, x2, x3;
	double y0, y1, y2, y3;
	double z0, z1, z2, z3;
	double t0, t1, t2, t3;
	double deltax, deltay, deltaz;
	double deltax1, deltay1, deltaz1;
	double deltax2, deltay2, deltaz2;
	double deltax3, deltay3, deltaz3;
	double Cx, Cy, Cz; // cross product in Lorentz force
	double gama;
	double RKx1, RKx2, RKx3, RKx4;
	double RKy1, RKy2, RKy3, RKy4;
	double RKz1, RKz2, RKz3, RKz4;
	double RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4;
	double RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4;
	double RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4;

    double B_tot;
//	double B_back; // background B-field
	
	x0 = *x; y0 = *y; z0 = *z;
	vx0 = *vx; vy0 = *vy; vz0 = *vz;
	t0 = *t;
	getdelta(&deltax, &deltay, &deltaz, vx0, vy0, vz0);	

// ===== now RK1 ====
//
    getEM(x0, y0, z0, t0, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);
    //B_tot = sqrt(B_totx * B_totx + B_toty * B_toty + B_totz * B_totz);
    /*
    if(iBessel == 0) {
        dt = epsilon / B_tot; // change dt for each position dt = epsilon * m / (q*B)
    }
    */
    h = dt / 2.0;
//    printf("%19.18e\n", dt);
//    printf("%10.9e\n", E_totz);

	cross_product(vx0, vy0, vz0, B_totx, B_toty, B_totz,
			&Cx, &Cy, &Cz);
	RKx1 = vx0 * c0 / L0;
	RKy1 = vy0 * c0 / L0;
	RKz1 = vz0 * c0 / L0;
//    printf("E_totx, E_toty, E_totz %15.14e %15.14e %15.14e\n", E_totx, E_toty, E_totz);
	RKdeltax1 = ChargeMass*(E_totx + Cx);
	RKdeltay1 = ChargeMass*(E_toty + Cy);
	RKdeltaz1 = ChargeMass*(E_totz + Cz);

//    printf("RK1 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);

// ==== now RK2 ====

	x1 = x0 + h * RKx1;
	y1 = y0 + h * RKy1;
	z1 = z0 + h * RKz1;
	t1 = t0 + h;

	deltax1 = deltax + h * RKdeltax1;
	deltay1 = deltay + h * RKdeltay1;
	deltaz1 = deltaz + h * RKdeltaz1;
	
	gama = getgamma(deltax1, deltay1, deltaz1);
	vx1 = deltax1 / gama;
	vy1 = deltay1 / gama;
	vz1 = deltaz1 / gama;

    getEM(x1, y1, z1, t1, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);

	cross_product(vx1, vy1, vz1, B_totx, B_toty, B_totz,
			&Cx, &Cy, &Cz);
	RKx2 = vx1 * c0 / L0;
	RKy2 = vy1 * c0 / L0;
	RKz2 = vz1 * c0 / L0;
	RKdeltax2 = ChargeMass*(E_totx + Cx);
	RKdeltay2 = ChargeMass*(E_toty + Cy);
	RKdeltaz2 = ChargeMass*(E_totz + Cz);

//    printf("RK2 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);
// ==== now RK3 ====

	x2 = x0 + h*RKx2;
	y2 = y0 + h*RKy2;
	z2 = z0 + h*RKz2;
	t2 = t0 + h;

	deltax2 = deltax + h*RKdeltax2;
	deltay2 = deltay + h*RKdeltay2;
	deltaz2 = deltaz + h*RKdeltaz2;
	
	gama = getgamma(deltax2, deltay2, deltaz2);
	vx2 = deltax2 / gama;
	vy2 = deltay2 / gama;
	vz2 = deltaz2 / gama;

    getEM(x2, y2, z2, t2, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);

	cross_product(vx2, vy2, vz2, B_totx, B_toty, B_totz,
			&Cx, &Cy, &Cz);
	RKx3 = vx2 * c0 / L0;
	RKy3 = vy2 * c0 / L0;
	RKz3 = vz2 * c0 / L0;
	RKdeltax3 = ChargeMass*(E_totx +Cx);
	RKdeltay3 = ChargeMass*(E_toty +Cy);
	RKdeltaz3 = ChargeMass*(E_totz +Cz);

//    printf("RK3 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);
// ==== now RK4 ====

	x3 = x0 + 2.0*h*RKx3;
	y3 = y0 + 2.0*h*RKy3;
	z3 = z0 + 2.0*h*RKz3;
	t3 = t0 + 2.0*h;

	deltax3 = deltax + 2.0*h*RKdeltax3;
	deltay3 = deltay + 2.0*h*RKdeltay3;
	deltaz3 = deltaz + 2.0*h*RKdeltaz3;
	
	gama = getgamma(deltax3, deltay3, deltaz3);
	vx3 = deltax3 / gama;
	vy3 = deltay3 / gama;
	vz3 = deltaz3 / gama;

    getEM(x3, y3, z3, t3, iBessel, &B_totx, &B_toty, &B_totz, &E_totx, &E_toty, &E_totz);

	cross_product(vx3, vy3, vz3, B_totx, B_toty, B_totz,
			&Cx, &Cy, &Cz);
	RKx4 = vx3 * c0 / L0;
	RKy4 = vy3 * c0 / L0;
	RKz4 = vz3 * c0 / L0;
	RKdeltax4 = ChargeMass*(E_totx +Cx);
	RKdeltay4 = ChargeMass*(E_toty +Cy);
	RKdeltaz4 = ChargeMass*(E_totz +Cz);

//    printf("RK4 Cx, Cy, Cz %10.8e %10.8e %10.8e\n ", Cx, Cy, Cz);

//==== finishing RK method. old (x,y,z,vx,vy,vz,t) -> new (x,y,z,vx,vy,vz,t)
/*        
        printf("RKx1, RKx2, RKx3, RKx4: %10.9e %10.9e %10.9e %10.9e\n",
		RKx1, RKx2, RKx3, RKx4);
        printf("RKy1, RKy2, RKy3, RKy4: %10.9e %10.9e %10.9e %10.9e\n",
		RKy1, RKy2, RKy3, RKy4);
        printf("RKz1, RKz2, RKz3, RKz4: %10.9e %10.9e %10.9e %10.9e\n",
		RKz1, RKz2, RKz3, RKz4);
        printf("RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4: %10.9e %10.9e %10.9e %10.9e\n",
		RKdeltax1, RKdeltax2, RKdeltax3, RKdeltax4);
        printf("RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4: %10.9e %10.9e %10.9e %10.9e\n",
		RKdeltay1, RKdeltay2, RKdeltay3, RKdeltay4);
        printf("RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4: %10.9e %10.9e %10.9e %10.9e\n",
		RKdeltaz1, RKdeltaz2, RKdeltaz3, RKdeltaz4);
*/
	*x += h/3.0 * ( RKx1 + 2.0*RKx2 + 2.0*RKx3 + RKx4);
	*y += h/3.0 * ( RKy1 + 2.0*RKy2 + 2.0*RKy3 + RKy4);
	*z += h/3.0 * ( RKz1 + 2.0*RKz2 + 2.0*RKz3 + RKz4);

	deltax += h/3.0 * (RKdeltax1 + 2.0*RKdeltax2 + 2.0*RKdeltax3 + RKdeltax4);
	deltay += h/3.0 * (RKdeltay1 + 2.0*RKdeltay2 + 2.0*RKdeltay3 + RKdeltay4);
	deltaz += h/3.0 * (RKdeltaz1 + 2.0*RKdeltaz2 + 2.0*RKdeltaz3 + RKdeltaz4);

	gama = getgamma(deltax, deltay, deltaz);

	*vx = deltax / gama;
	*vy = deltay / gama;
	*vz = deltaz / gama;
	*t += dt;

}

/******************************************************************************
 * Calculate gama * beta. gama is Lorentz factor. Beta is the ratio of velocity
 * over speed of light.
 ******************************************************************************/
void getdelta(double *delta_x, double *delta_y, double *delta_z,
		double beta_x, double beta_y, double beta_z)
{
	double beta_mag, gama;
	beta_mag = sqrt(beta_x * beta_x + beta_y * beta_y + beta_z * beta_z);
	gama = 1.0/sqrt(1.0 - beta_mag * beta_mag);
//    printf("gama %20.19e\n", gama);
	*delta_x = gama * beta_x;
	*delta_y = gama * beta_y;
	*delta_z = gama * beta_z;
}

/******************************************************************************
 * Calculate gama from gama * beta.
 ******************************************************************************/
double getgamma(double dx, double dy, double dz)
{
	double delta, gamma;

	delta = sqrt(dx * dx + dy * dy + dz * dz);
    gamma = sqrt(1.0 + delta * delta);
//    printf("%20.19e\n ", gamma);
	return gamma;
}

/******************************************************************************
 * Swap functions.
 ******************************************************************************/
void swap_double(double *a, double *b)
{
    double temp;
    temp = *a;
    *a = *b;
    *b = temp;
}

void swap_int(int *a, int *b)
{
    int temp;
    temp = *a;
    *a = *b;
    *b = temp;
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

