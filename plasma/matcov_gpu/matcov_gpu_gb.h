#ifndef MATCOV_GPU_GB
#define MATCOV_GPU_GB


#include<cuda.h>
#include<cuda_runtime_api.h>

// tomo_gpu_struct -> gtomo_structu 
struct gtomo_struct {
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
  long max_Nl0;
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

void process_err(cudaError_t e, const char* str);
void matts_gpu_gb(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, struct tomo_struct tomo);
void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu);
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu,struct tomo_struct tomo);
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu,struct tomo_struct tomo);
void matcov_gpu_3(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
void matcov_gpu_4(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu);
#endif // MATCOV_GPU_GB
