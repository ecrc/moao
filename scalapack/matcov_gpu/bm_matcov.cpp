#include <sys/time.h>
#include<stdio.h>
#include<cuda.h>
#include<cuda_runtime_api.h>
#include<cublas.h>
#include "matcov.h"
#include "matcov_gpu_gb.h"
//#include <cfitsio/fitsio2.h>

void write_sys_file(int nssp)
{
	FILE *fp = fopen("sys-params.txt", "w");
	if(fp == NULL){printf("error opening file\n"); exit(1);}
	
	fprintf(fp, "tel diam in meters\n");
	fprintf(fp, "40.0\n");
	fprintf(fp, "cent obs\n");
	fprintf(fp, "0.2\n");
	fprintf(fp, "N WFS\n");
	fprintf(fp, "10\n");
	fprintf(fp, "Nssp n subaps per wfs\n");
	fprintf(fp, "%d\n", nssp);
	fprintf(fp, "GsAlt\n");
	fprintf(fp, "1e-05 1e-05 1e-05 1e-05 1e-05 1e-05 0 0 0 0\n");
	fprintf(fp, "type\n");
	fprintf(fp, "2 2 2 2 2 2 1 1 1 1\n");
	fprintf(fp, "nlgs\n");
	fprintf(fp, "6\n");
	fprintf(fp, "alphaX in arcsec\n");
	fprintf(fp, "102.0 51.0 -51.0 -102.0 -51.0 51.0 -16.0 -32.0 52.0 0.0\n");
	fprintf(fp, "alphaY in arcsec\n");
	fprintf(fp, "0.0 88.3346 88.3346 0.0 -88.3346 -88.3346 -36.0 38.0 -5.0\n");
	fprintf(fp, "XPup\n");
	fprintf(fp, "0 0 0 0 0 0 0 0 0 0\n");
	fprintf(fp, "YPup\n");
	fprintf(fp, "0 0 0 0 0 0 0 0 0 0\n");
	fprintf(fp, "thetaML\n");
	fprintf(fp, "0 0 0 0 0 0 0 0 0 0\n");
	fprintf(fp, "thetaCam\n");
	fprintf(fp, "0 0 0 0 0 0 0 0 0 0\n");
	fprintf(fp, "sensibilite\n");
	fprintf(fp, "1 1 1 1 1 1 1 1 1 1\n");
	fprintf(fp, "tracking\n");
	fprintf(fp, "0 0 0\n");
	fprintf(fp, "pasDphi\n");
	fprintf(fp, "0.0001\n");
	fprintf(fp, "ncpu\n");
	fprintf(fp, "12\n");
	fprintf(fp, "lgs_cst\n");
	fprintf(fp, "0\n");
	fprintf(fp, "noise_var\n");
	fprintf(fp, "0.0\n");
	fprintf(fp, "spot_width arcsec\n");
	fprintf(fp, "1\n");
	fprintf(fp, "lgs_alt (m)\n");
	fprintf(fp, "90000.0\n");
	fprintf(fp, "lgs_depth (m)\n");
	fprintf(fp, "10000.0\n");

	fclose(fp);
}

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
	if(argc < 4)
	{
		printf("USAGE: %s <nssp-start> <nssp-end> <nssp-step>\n", argv[0]);
		exit(1);
	}
	
	int nssp_start = atoi(argv[1]);
	int nssp_end = atoi(argv[2]);
	int nssp_step = atoi(argv[3]);
	
	int type_mat[3] = {1, 3, 4};
	
	struct tomo_struct tomo;
	
	for(int t = 0; t < 3; t++)
	{
		int type = type_mat[t];
		
		printf("Matrix Name       M           N       cpu time (ms)   gpu time (ms)   speedup   max error \n");
  		printf("-----------   ---------   ---------   -------------   -------------   -------   --------- \n");
  			
		for(int n = nssp_start; n <= nssp_end; n += nssp_step)
		{
			write_sys_file(n);
			
			//initialize tomo struct
			init_tomo_sys(&tomo);
  			init_tomo_atm(&tomo);
  			
  			long nsubaps_offaxis;       // number of sub-aperture of offaxis WFSs
  			nsubaps_offaxis = 0;
  			for (int i = 0; i < tomo.Nw-1; i++)
    			nsubaps_offaxis += tomo.Nsubap[i];
			
			
  			int dim, size_cpu, size_gpu;
  			dim = (nsubaps_offaxis * 2);
  			int M, N, LDA, LDA_;
  			
			struct timeval start_cpu, stop_cpu;
			cudaEvent_t start, stop;
			cudaEventCreate(&start);
  			cudaEventCreate(&stop);
  		
  			// init on gpu
  			struct gtomo_struct tomo_gpu;
  			init_tomo_gpu_gb(&tomo_gpu,tomo);
  			
  			update_tomo_sys_gpu_gb(&tomo_gpu,tomo); // done once
  			cudaStreamSynchronize(tomo_gpu.matcov_stream);
		
			// at least this has to run once
			// this is done everytime we recompute cmm
			update_tomo_atm_gpu_gb(&tomo_gpu,tomo); 
			cudaStreamSynchronize(tomo_gpu.matcov_stream);
			
  			double cpu_time;
  			float gpu_time;
  			double max_err;
  			int index = 0;
  			
  			int xoffset, yoffset;
  			cublasStatus_t s ;
  			
  			// Cmm generation
  			if(type == 1)
  			{
  				M = dim; 
  				N = dim;
				
  				LDA = M;
  				LDA_ = ((M+31)/32)*32;
				
  				size_cpu = N * LDA;
  				size_gpu = N * LDA_;
			
  				// host-side matrices
  				double* Cmm 	= (double*)malloc(size_cpu * sizeof(double));
  				double* Cmm_gpu	= (double*)malloc(size_cpu * sizeof(double));
				
				printf("Cmm           ");
  				printf("%-9d   ", M);
  				printf("%-9d   ", N);
  		
  				// cpu test
  				gettimeofday(&start_cpu, NULL);
  				fill_Caa(tomo, Cmm);
  				gettimeofday(&stop_cpu, NULL);
				
				cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;
				printf("%-13.2f   ", cpu_time);
			
  				// gpu test
  				xoffset = 0;
  				yoffset = 0;
  				
  				double *da;
  				cudaMalloc((void**)&da, size_gpu * sizeof(double));
  				if(da == NULL){printf("error allocating da on gpu\n");exit(1);}
					
  				cudaEventRecord(start, 0);
  				
  				// compute Cmm
  				update_tomo_atm_gpu_gb(&tomo_gpu,tomo); // this is done everytime we recompute cmm
  				cudaStreamSynchronize(tomo_gpu.matcov_stream); 
  				tomo.part = 1;
  				matcov_gpu_4(da, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);
  				cudaStreamSynchronize(tomo_gpu.matcov_stream);
  				//---
  		
  				cudaEventRecord(stop, 0);
  				cudaEventSynchronize(stop);
  				cudaEventElapsedTime(&gpu_time, start, stop);
		
				printf("%-13.2f   ", gpu_time);
  				printf("%-7.2f   ", (cpu_time/gpu_time));
  				 
  				s = cublasGetMatrix(M, N, sizeof(double), da, LDA_, Cmm_gpu, LDA);
  				if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
  				
  				// test correctness
  				{
  					double* ref = Cmm;
  					double* res = Cmm_gpu;
  					max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  				}
  				
  				printf("%-9.7e  @ %d", max_err, index);
  				printf("\n");
  				
  				// free resources
  				if(Cmm)free(Cmm);
  				if(Cmm_gpu)free(Cmm_gpu);
  				if(da)cudaFree(da);
  			} // end of Cmm generation test
  			
  			
  			// Cpm test
  			if(type == 3)
  			{
  				M = tomo.Nsubap[tomo.Nw-1] * 2; 
  				N = dim;
  				
  				LDA = M;
  				LDA_ = ((M+31)/32)*32;
				
				size_cpu = N * LDA;
				size_gpu = N * LDA_;
		
				// host-side matrices
  				double* Cpm 	= (double*)malloc(size_cpu * sizeof(double));
  				double* Cpm_gpu	= (double*)malloc(size_cpu * sizeof(double));
		
				printf("Cpm           ");
				printf("%-9d   ", M);
				printf("%-9d   ", N);
  		
  				// cpu test
  				gettimeofday(&start_cpu, NULL);
  				fill_Cmaa(tomo, Cpm);
  				gettimeofday(&stop_cpu, NULL);
  				
  				cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;
				
				printf("%-13.2f   ", cpu_time);
		
				double *db;
  				cudaMalloc((void**)&db, size_gpu * sizeof(double));
  				if(db == NULL){printf("error allocating db on gpu\n");exit(1);}
		
  				cudaEventRecord(start, 0);
  				
  				tomo.part = 3;  
  				matcov_gpu_4(db, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);
  				cudaStreamSynchronize(tomo_gpu.matcov_stream);
  		
  				cudaEventRecord(stop, 0);
  				cudaEventSynchronize(stop);
  				cudaEventElapsedTime(&gpu_time, start, stop);
  		
  				printf("%-13.2f   ", gpu_time);
  				printf("%-7.2f   ", (cpu_time/gpu_time));
  				
  				s = cublasGetMatrix(M, N, sizeof(double), db, LDA_, Cpm_gpu, LDA);
  				if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
  				// test correctness
  				{
    				double* ref = Cpm;
    				double* res = Cpm_gpu;
    				max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  				}
  				
  				printf("%-9.7e  @ %d", max_err, index);
				printf("\n");
				
				// free resources
				if(Cpm)free(Cpm);
  				if(Cpm_gpu)free(Cpm_gpu);
  				if(db)cudaFree(db);
  			}	
			
			// Cpp test
			if(type == 4)
			{
				M = tomo.Nsubap[tomo.Nw-1] * 2; 
				N = tomo.Nsubap[tomo.Nw-1] * 2;
			
				LDA = M;
				LDA_ = ((M+31)/32)*32;
			
  				size_cpu = N * LDA;
  				size_gpu = N * LDA_;
  				
  				// host-side matrices
  				double* Cpp 	= (double*)malloc(size_cpu * sizeof(double));
  				double* Cpp_gpu	= (double*)malloc(size_cpu * sizeof(double));
			
  				printf("Cpp           ");
  				printf("%-9d   ", M);
  				printf("%-9d   ", N);
  		
  				// cpu test
  				gettimeofday(&start_cpu, NULL);
  				matcov_cpp(tomo, Cpp);
  				gettimeofday(&stop_cpu, NULL);
			
  				cpu_time = (stop_cpu.tv_sec - start_cpu.tv_sec) * 1000.0 + (stop_cpu.tv_usec - start_cpu.tv_usec) * 0.001;
		
  				printf("%-13.2f   ", cpu_time);
  		
  				double *dc;
  				cudaMalloc((void**)&dc, size_gpu * sizeof(double));
  				if(dc == NULL){printf("error allocating dc on gpu\n");exit(1);}
  				
  				cudaEventRecord(start, 0);
				
				tomo.part = 4;
				matts_gpu_gb(dc, M, N, xoffset, yoffset, LDA_, tomo, &tomo_gpu);
				cudaStreamSynchronize(tomo_gpu.matcov_stream);
				
				cudaEventRecord(stop, 0);
				cudaEventSynchronize(stop);
				cudaEventElapsedTime(&gpu_time, start, stop);
		
  				printf("%-13.2f   ", gpu_time);
  				printf("%-7.2f   ", (cpu_time/gpu_time));
  			
  				s = cublasGetMatrix(M, N, sizeof(double), dc, LDA_, Cpp_gpu, LDA);
  				if(s != CUBLAS_STATUS_SUCCESS)printf("error %d in cublasGetMatrix \n", (int)s);
  				
  				// test correctness
  				{
    				double* ref = Cpp;
    				double* res = Cpp_gpu;
    				max_err = compute_max_error(ref, res, size_cpu, 1, &index);
  				}
  				printf("%-9.7e  @ %d", max_err, index);
				printf("\n");
				
				// free resources
				if(dc)cudaFree(dc);
  				if(Cpp)free(Cpp);
  				if(Cpp_gpu)free(Cpp_gpu);
			}
			
			// free tomo
			free_tomo_gpu_gb(&tomo_gpu);
			free_tomo(&tomo);
			
			cudaEventDestroy(start);
			cudaEventDestroy(stop);
		}
		
		printf("\n\n");
	}
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
