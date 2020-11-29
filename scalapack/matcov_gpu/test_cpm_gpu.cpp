#include <sys/time.h>
#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime_api.h>
#include<cublas.h>
#include "matcov.h"
#include "matcov_gpu.h"

double compute_max_error(double* ref, double *res, int n, int inc, int* index)
{
	int i;
	double max_err = -1.0;
	double err = -1.0;
	inc = abs(inc);
	int ind = 0;
	for(i = 0; i < n; i++)
	{
		// ignore values lower than epsilon
		if(fabs(ref[i*inc]) < 1e-16)continue;
		
		err = fabs(res[i * inc] - ref[i * inc]);
		if(ref[i * inc] != 0.0)err /= fabs(ref[i * inc]);
		if(err > max_err)
		{
			max_err = err;
			ind = i;
		}
		//if (i > 100 && i < 200)
		//{
		//	printf("[%-4d]: %+-.4e   %+-.4e   %e", i, ref[i * inc], res[i * inc], err);
		//	if(err > 1e-12)printf("    ---> error");
		//	printf("\n");
		//}
		
	}
	*index = ind;
	return max_err;
}

//---------------------------------------------------------------------------
int main(int argc, char* argv[])
{
	if(argc < 2)
	{
		printf("USAGE: %s <M> \n", argv[0]);
		printf("====================================================\n");
		printf("M : nssp\n");
		return -1;
	}
	
	
	int nssp = atoi(argv[1]);

	struct tomo_struct tomo;
  	initialize_tomo2(&tomo,nssp);
  
  	//print_tomo(tomo);

	long nsubaps_offaxis;       // number of sub-aperture of offaxis WFSs
	nsubaps_offaxis = 0;
	for (int i = 0; i < tomo.Nw-1; i++)
	  nsubaps_offaxis += tomo.Nsubap[i];


	int dim, size_cpu, size_gpu;
	dim = (nsubaps_offaxis * 2);
	int N = 2*tomo.Nsubap[tomo.Nw-1]; 
	int M = dim;
	
	int LDA = M;
	int LDA_ = ((M+31)/32)*32;
	
	size_cpu = N * LDA;
	size_gpu = N * LDA_;
	
	// host-side matrices
	double* Caa 	= (double*)malloc(size_cpu * sizeof(double));
	double* Caa_gpu	= (double*)malloc(size_cpu * sizeof(double));
	
	printf("M = %d, N = %d, LDA = %d, LDA_ = %d \n", M, N, LDA, LDA_);
	printf("cpu test \n");
	struct timeval start_cpu, stop_cpu;
	gettimeofday(&start_cpu, NULL);
	fill_Cmaa_dam(tomo, Caa);
	gettimeofday(&stop_cpu, NULL);
	
	double cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;

	int xoffset = 0;
	int yoffset = 0;
	int print = 0;
	
	printf("gpu test \n");
	double *da;
	cudaMalloc((void**)&da, size_gpu * sizeof(double));
	if(da == NULL){printf("error allocating da on gpu\n");exit(1);}
	
	cudaEvent_t start, stop;
	cudaEventCreate(&start);
	cudaEventCreate(&stop);
	
	cudaEventRecord(start, 0);
	tomo.part = 3;		// fill Cmaa matrix
	matcov_gpu(da, M, N, xoffset, yoffset, LDA_, tomo);
	cudaEventRecord(stop, 0);
	cudaEventSynchronize(stop);
	
	float time;
	cudaEventElapsedTime(&time, start, stop);
	
	cublasStatus_t s = cublasGetMatrix(M, N, sizeof(double), da, LDA_, Caa_gpu, LDA);
	if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
	
	printf("\n======================\n");
	printf("GPU generation time = %.2f ms", time);
	printf("\n======================\n");
	printf("\n======================\n");
	printf("CPU time = %.2f ms", cpu_time);
	printf("\n======================\n");
	
	double max_err;
	int index = 0;
	// test correctness
	{
		double* ref = Caa;
		double* res = Caa_gpu;
		max_err = compute_max_error(ref, res, size_cpu, 1, &index);
	}
	printf("max error = %e @ [%d]\n", max_err, index);
	printf("======================\n");

	// free resources
	if(Caa)free(Caa);
	if(Caa_gpu)free(Caa_gpu);
	if(da)cudaFree(da);
	
	cudaEventDestroy(start);
	cudaEventDestroy(stop);
}
	

