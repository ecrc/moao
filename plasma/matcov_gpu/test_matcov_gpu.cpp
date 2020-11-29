#include <sys/time.h>
#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime_api.h>
#include<cublas.h>
#include "matcov.h"
#include "matcov_gpu.h"
//#include <cfitsio/fitsio2.h>

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
  /*
    if(argc < 2)
    {
    printf("USAGE: %s <M> \n", argv[0]);
    printf("====================================================\n");
    printf("M : nssp\n");
    return -1;
    }
	
	
    int nssp = atoi(argv[1]);
  */
  struct tomo_struct tomo;
  //initialize_tomo2(&tomo,nssp);
  init_tomo_sys(&tomo);
  init_tomo_atm(&tomo);
  
  //print_tomo(tomo);

  long nsubaps_offaxis;       // number of sub-aperture of offaxis WFSs
  nsubaps_offaxis = 0;
  for (int i = 0; i < tomo.Nw-1; i++)
    nsubaps_offaxis += tomo.Nsubap[i];


  int dim, size_cpu, size_gpu;
  dim = (nsubaps_offaxis * 2);
  int M = dim; 
  int N = dim;
	
  int LDA = M;
  int LDA_ = ((M+31)/32)*32;
	
  size_cpu = N * LDA;
  size_gpu = N * LDA_;
	
  // host-side matrices
  double* Cmm 	= (double*)malloc(size_cpu * sizeof(double));
  double* Cmm_gpu	= (double*)malloc(size_cpu * sizeof(double));
	
  printf("M = %d, N = %d, LDA = %d, LDA_ = %d \n", M, N, LDA, LDA_);
  printf("cpu test \n");
  struct timeval start_cpu, stop_cpu;
  gettimeofday(&start_cpu, NULL);
  fill_Caa(tomo, Cmm);
  gettimeofday(&stop_cpu, NULL);
	
  double cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;

  printf("\n======================\n");
  printf("CPU time = %.2f ms", cpu_time);
  printf("\n======================\n");
	
  printf("gpu test");
  int xoffset = 0;
  int yoffset = 0;
  int print = 0;
	
  //  	initialize_tomo(&tomo);
  //print_tomo(tomo);

  double *da;
  cudaMalloc((void**)&da, size_gpu * sizeof(double));
  if(da == NULL){printf("error allocating da on gpu\n");exit(1);}
	
  cudaEvent_t start, stop;
  float time;
  /*
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
	
    cudaEventRecord(start, 0);
    tomo.part = 1;		// fill Cmm matrix
    matcov_gpu(da, M, N, xoffset, yoffset, LDA_, tomo);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
	
    cudaEventElapsedTime(&time, start, stop);
	
    printf("\n======================\n");
    printf("GPU1 time = %.2f ms", time);
    printf("\n======================\n");
  */
  struct tomo_gpu_struct tomo_gpu;
  init_tomo_gpu(&tomo_gpu,tomo);

  update_tomo_sys(&tomo_gpu,tomo); // this is done once

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
	
  cudaEventRecord(start, 0);

  update_tomo_atm(&tomo_gpu,tomo); // this is done everytime we recompute cmm

  tomo.part = 1;

  matcov_gpu4(da, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);

  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  printf("\n======================\n");
  printf("GPU2 time = %.2f ms", time);
  printf("\n======================\n");
	
  cublasStatus_t s = cublasGetMatrix(M, N, sizeof(double), da, LDA_, Cmm_gpu, LDA);
  if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
	
  double max_err;
  int index = 0;
  // test correctness
  {
    double* ref = Cmm;
    double* res = Cmm_gpu;
    max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  }
  printf("max error = %e @ [%d]\n", max_err, index);
  printf("======================\n");


  
  M = tomo.Nsubap[tomo.Nw-1] * 2; 
  N = dim;
	
  LDA = M;
  LDA_ = ((M+31)/32)*32;
	
  size_cpu = N * LDA;
  size_gpu = N * LDA_;

  // host-side matrices
  double* Cpm 	= (double*)malloc(size_cpu * sizeof(double));
  double* Cpm_gpu	= (double*)malloc(size_cpu * sizeof(double));
	
  printf("Cpm cpu test \n");
  printf("M = %d, N = %d, LDA = %d, LDA_ = %d \n", M, N, LDA, LDA_);
  gettimeofday(&start_cpu, NULL);
  fill_Cmaa(tomo, Cpm);
  gettimeofday(&stop_cpu, NULL);
	
  cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;

  printf("\n======================\n");
  printf("CPU time = %.2f ms", cpu_time);
  printf("\n======================\n");

  double *db;
  cudaMalloc((void**)&db, size_gpu * sizeof(double));
  if(db == NULL){printf("error allocating db on gpu\n");exit(1);}
  tomo.part = 3;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
	
  cudaEventRecord(start, 0);

  //update_tomo_atm(&tomo_gpu,tomo); // this is done everytime we recompute cmm

  matcov_gpu4(db, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);
  //matcov_gpu(db, M, N, xoffset, yoffset, LDA_, tomo);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  printf("\n======================\n");
  printf("GPU2 time = %.2f ms", time);
  printf("\n======================\n");

  s = cublasGetMatrix(M, N, sizeof(double), db, LDA_, Cpm_gpu, LDA);
  if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
  // test correctness
  {
    double* ref = Cpm;
    double* res = Cpm_gpu;
    max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  }
  printf("max error = %e @ [%d]\n", max_err, index);
  printf("======================\n");


  M = tomo.Nsubap[tomo.Nw-1] * 2; 
  N = tomo.Nsubap[tomo.Nw-1] * 2;
	
  LDA = M;
  LDA_ = ((M+31)/32)*32;
	
  size_cpu = N * LDA;
  size_gpu = N * LDA_;

  // host-side matrices
  double* Cpp 	= (double*)malloc(size_cpu * sizeof(double));
  double* Cpp_gpu	= (double*)malloc(size_cpu * sizeof(double));
	
  printf("Cpp cpu test \n");
  printf("M = %d, N = %d, LDA = %d, LDA_ = %d \n", M, N, LDA, LDA_);
  gettimeofday(&start_cpu, NULL);
  matcov_cpp(tomo, Cpp);
  gettimeofday(&stop_cpu, NULL);
	
  cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;

  printf("\n======================\n");
  printf("CPU time = %.2f ms", cpu_time);
  printf("\n======================\n");
  double *dc;
  cudaMalloc((void**)&dc, size_gpu * sizeof(double));
  if(dc == NULL){printf("error allocating dc on gpu\n");exit(1);}
  tomo.part = 4;

  cudaEventCreate(&start);
  cudaEventCreate(&stop);
	
  cudaEventRecord(start, 0);

  //update_tomo_atm(&tomo_gpu,tomo); // this is done everytime we recompute cmm

  matts_gpu(dc, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);
  //matcov_gpu(db, M, N, xoffset, yoffset, LDA_, tomo);
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&time, start, stop);

  printf("\n======================\n");
  printf("GPU2 time = %.2f ms", time);
  printf("\n======================\n");

  s = cublasGetMatrix(M, N, sizeof(double), dc, LDA_, Cpp_gpu, LDA);
  if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
  // test correctness
  {
    double* ref = Cpp;
    double* res = Cpp_gpu;
    max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  }
  printf("max error = %e @ [%d]\n", max_err, index);
  printf("======================\n");


 // free resources
  if(dc)cudaFree(dc);
  if(Cpp)free(Cpp);
  if(Cpp_gpu)free(Cpp_gpu);

  if(Cpm)free(Cpm);
  if(Cpm_gpu)free(Cpm_gpu);
  if(db)cudaFree(db);
  

  free_tomo_gpu(&tomo_gpu);

  if(Cmm)free(Cmm);
  if(Cmm_gpu)free(Cmm_gpu);
  if(da)cudaFree(da);
	
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
}
	

  /*
fitsfile *fptr; // FITS file pointer, defined in fitsio.h 
	int status = 0; // CFITSIO status value MUST be initialized to zero! 
fits_create_file(&fptr, "ref.fits", &status);
int naxis = 2;
//long naxis_ipct[] = { nsubaps_offaxis * 2, nsubaps_offaxis * 2};
 long naxis_ipct[] = {M  , N};
//long nelements = nsubaps_offaxis * 2 * nsubaps_offaxis * 2;
long nelements = M * N;
fits_create_img(fptr, DOUBLE_IMG, naxis, naxis_ipct, &status);
long fpixel[] = { 1, 1 };  // coordinate in each dimension of the first pixel to be written 
fits_write_pix(fptr, TDOUBLE, fpixel, nelements, Cpm, &status);
fits_close_file(fptr, &status);
fits_report_error(stderr, status);
	
fitsfile *fptr2; 
fits_create_file(&fptr2, "res.fits", &status);
fits_create_img(fptr2, DOUBLE_IMG, naxis, naxis_ipct, &status);
fits_write_pix(fptr2, TDOUBLE, fpixel, nelements, Caa_dam, &status);
fits_close_file(fptr2, &status);
fits_report_error(stderr, status);
  */
