#include <sys/time.h>
#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime_api.h>
#include<cublas.h>
#include "matcov.h"
#include "matcov_gpu.h"

double compute_max_error(double* ref, double *res, int n, int inc)
{
	int i;
	double max_err = -1.0;
	double err = -1.0;
	inc = abs(inc);
	for(i = 0; i < n; i++)
	{
		err = fabs(res[i * inc] - ref[i * inc]);
		if(ref[i * inc] != 0.0)err /= fabs(ref[i * inc]);
		if(err > max_err)max_err = err;
	}
	return max_err;
}
//------------------------------------------------------------------------------------

int main(int argc, char* argv[])
{
	if(argc < 3)
	{
		printf("USAGE: %s <device-id> <nssp> \n", argv[0]);
		return -1;
	}
	
	int gpu_id = atoi(argv[1]);
	int nssp = atoi(argv[2]);
	
	cudaError_t e;
	
	e = cudaSetDevice(gpu_id);
	process_error(e, "set device");
	
	struct tomo_struct tomo;
  	initialize_tomo2(&tomo, nssp);
  
  	print_tomo(tomo);
  	
	// %%%%%%% Pre-computation of DPHI %%%%%%%%%%
	//Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
	const long cNw = tomo.Nw;
	const long cNlayer = tomo.Nlayer;
	const double crmax = tomo.rmax;
	const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
	const long cNdphi = floor(crmax*pasDPHI)+1;
	const double convert = (double)(cNdphi-1)/(crmax+1./pasDPHI);
	double rr[cNdphi];
	long i; 
	
	long indexL0[cNlayer]; //link between index in L0 and index in L0diff
	double **tabDPHI;
	double *tabDPHI_d;
	
	//rr varie de 0 à rmax+pasDphi (évite de traiter à part les cas où rr =rmax dans la fonction de DPHI)
	for (i=0;i<cNdphi;i++) rr[i] = (double)(i)/ convert;
	
	double* rr_d;
	cudaMalloc((void**)&rr_d, cNdphi*sizeof(double));
	
	printf("cpu test \n");
	// cpu generation
	struct timeval start_cpu, stop_cpu;
	gettimeofday(&start_cpu, NULL);
	tabDPHI = tabulateDPHI(rr, tomo, cNdphi, indexL0);
	gettimeofday(&stop_cpu, NULL);
	
	double cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;
	
	printf("gpu test \n");
	// gpu generation
	cudaEvent_t start_gpu, stop_gpu;
	cudaEventCreate(&start_gpu);
	cudaEventCreate(&stop_gpu);
	
	int Nl0_; // used to know the size of the array
	
	cudaEventRecord(start_gpu, 0);
	cudaMemcpy(rr_d, rr, cNdphi*sizeof(double), cudaMemcpyHostToDevice);
	tabDPHI_d = tabulateDPHI_gpu(tomo, cNdphi, indexL0, (int*)&Nl0_,convert);
	cudaEventRecord(stop_gpu, 0);
	cudaEventSynchronize(stop_gpu);
	
	float gpu_time;
	cudaEventElapsedTime(&gpu_time, start_gpu, stop_gpu);
	
	double max_err;
	int size;
	// test correctness
	{
		size = cNdphi * Nl0_;
		double* tmp = (double*)malloc(size * sizeof(double));
		cudaMemcpy(tmp, tabDPHI_d, size * sizeof(double), cudaMemcpyDeviceToHost);
		double* ref = tabDPHI[0];
		double* res = tmp;
		max_err = compute_max_error(ref, res, size, 1);
	}
	
	printf("cNdphi = %lu, Nl0 = %d, total size = %f KB", cNdphi, Nl0_, size*sizeof(double)/1024.0);
	printf("\n======================\n");
	printf("cpu generation time = %.2f ms\n", cpu_time);
	printf("gpu generation time = %.2f ms\n", gpu_time);
	if( (cpu_time/gpu_time) >= 1)
		printf("GPU is faster by speedup %.2f\n", (cpu_time/gpu_time) );
	else
		printf("CPU is faster by speedup %.2f\n", (gpu_time/cpu_time) );
	
	printf("max error = %e \n", max_err);
	printf("======================\n");
	
	// free resources
	if(rr_d)cudaFree(rr_d);
	if(tabDPHI[0])free(tabDPHI[0]);
	if(tabDPHI)free(tabDPHI);
	if(tabDPHI_d)cudaFree(tabDPHI_d);

	// destroy events
	cudaEventDestroy(start_gpu);
	cudaEventDestroy(stop_gpu);
}
	

