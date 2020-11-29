
#ifndef MATCOV_GPU
#define MATCOV_GPU
// header for functions

#include<cuda.h>
#include<cuda_runtime_api.h>

struct tomo_gpu_struct {
  long   *ioff_d;
  long   *Nssp_d;

  double *alphaX_d;
  double *alphaY_d;
  double *GsAlt_d;
  double *diamPup_d;
  double *thetaML_d;
  double *X_d;
  double *Y_d;
  double *XPup_d;
  double *YPup_d;

  long    Nl0;
  long   *indexL0_d;
  double *L0diff_d;
  double *h_d;
  double *cn2_d;
  double *tabDPHI_d;
  double *u_d;
  double *v_d;
  double *sspSizeL_d;
  /*
  double *Cmm_d;
  double *Cpm_d;
  double *R_d;
  */

  cudaStream_t matcov_stream;
};

void process_error(cudaError_t e, const char* str);

void matcov_gpu(double* data, int nrows, int ncols, 
				int xoffset, int yoffset, int lda , struct tomo_struct tomo);

void matcov_gpu2(double* data, int nrows, int ncols, 
				int xoffset, int yoffset, int lda , struct tomo_struct tomo);
				
void matcov_gpu3(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, 
		 struct tomo_gpu_struct *tomo_gpu);
void matts_gpu(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, 
		 struct tomo_gpu_struct *tomo_gpu);
void matcov_gpu4(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, 
		 struct tomo_gpu_struct *tomo_gpu);

double* tabulateDPHI_gpu(struct tomo_struct tomo, long Ndphi, long *indexL0, int* Nl0_, double convert);

void init_tomo_gpu(struct tomo_gpu_struct *tomo_gpu, struct tomo_struct tomo);
void free_tomo_gpu(struct tomo_gpu_struct *tomo_gpu);
void update_tomo_sys(struct tomo_gpu_struct *tomo_gpu,struct tomo_struct tomo);
void update_tomo_atm(struct tomo_gpu_struct *tomo_gpu,struct tomo_struct tomo);
void sub_pos_gpu(struct tomo_gpu_struct *tomo_gpu, struct tomo_struct tomo);
void tab_dphi_gpu(double *tab_dphi, struct tomo_struct tomo, struct tomo_gpu_struct *tomo_gpu, long Ndphi, double *L0diff_d, int Nl0, double convert);

#endif // MATCOV_GPU
