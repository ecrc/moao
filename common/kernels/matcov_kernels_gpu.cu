/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "matcov_kernels_gpu.hpp"

/*  Tuning parameters of tbulateDPHI kernel*/
#define tabDPHI_thread_x    (256)

/*  Tuning parameters of matcov GPU Kernel */
// Thread block size (x, y), 
// max #threads per block is 512 for fermi and 1024 for kepler
#define matcov_thread_x (8)
#define matcov_thread_y (8)

//============================================================================================
//============================= tabDPHI FUNCTIONS/KERNEL(s) ==================================
//============================================================================================
 template<typename T>
inline __device__ double macdo_x56_gpu(T x, int k)
/* DOCUMENT  macdo_x56_gpu(x)

 Computation of the function
 f(x) = x^(5/6)*K_{5/6}(x)
 using a series for the esimation of K_{5/6}, taken from Rod Conan thesis :
 K_a(x)=1/2 \sum_{n=0}^\infty \frac{(-1)^n}{n!}
 \left(\Gamma(-n-a) (x/2)^{2n+a} + \Gamma(-n+a) (x/2)^{2n-a} \right) ,
 with a = 5/6.

 Setting x22 = (x/2)^2, setting uda = (1/2)^a, and multiplying by x^a,
 this becomes :
 x^a * Ka(x) = 0.5 $ -1^n / n! [ G(-n-a).uda x22^(n+a) + G(-n+a)/uda x22^n ]
 Then we use the following recurrence formulae on the following quantities :
 G(-(n+1)-a) = G(-n-a) / -a-n-1
 G(-(n+1)+a) = G(-n+a) /  a-n-1
 (n+1)! = n! * (n+1)
 x22^(n+1) = x22^n * x22
 and at each iteration on n, one will use the values already computed at step (n-1).
 The values of G(a) and G(-a) are hardcoded instead of being computed.

 The first term of the series has also been skipped, as it
 vanishes with another term in the expression of Dphi.

 SEE ALSO:
 */
{
	const double a = 5. / 6.;
	const double x2a = pow(x, 2. * a), x22 = x * x / 4.;
	T x2n;               // x^2.a, etc
	T s = 0.0;
	int n;

	const double Ga[11] = { 0, 12.067619015983075, 5.17183672113560444,
    						0.795667187867016068, 0.0628158306210802181, 0.00301515986981185091,
    						9.72632216068338833e-05, 2.25320204494595251e-06, 3.93000356676612095e-08,
    						5.34694362825451923e-10, 5.83302941264329804e-12 };

	const double Gma[11] = { -3.74878707653729304, -2.04479295083852408,
    						-0.360845814853857083, -0.0313778969438136685, -0.001622994669507603,
    						-5.56455315259749673e-05, -1.35720808599938951e-06,
   							-2.47515152461894642e-08, -3.50257291219662472e-10,
    						-3.95770950530691961e-12, -3.65327031259100284e-14 };

	x2n = 0.5;                           // init (1/2) * x^0

	s = Gma[0] * x2a;
	s *= x2n;

	// prepare recurrence iteration for next step
	x2n *= x22;    // x^n

	#pragma unroll
	for (n = 1; n <= 10; n++)
	{
  		s += (Gma[n] * x2a + Ga[n]) * x2n;
  		// prepare recurrence iteration for next step
  		x2n *= x22;    // x^n
	}
	return s;
}
//------------------------------------------------------------------------------------
template<typename T>
inline __device__ double asymp_macdo_gpu(T x)
/* DOCUMENT asymp_macdo_gpu(x)

 Computes a term involved in the computation of the phase struct
 function with a finite outer scale according to the Von-Karman
 model. The term involves the MacDonald function (modified bessel
 function of second kind) K_{5/6}(x), and the algorithm uses the
 asymptotic form for x ~ infinity.
 Warnings :
 - This function makes a doubleing point interrupt for x=0
 and should not be used in this case.
 - Works only for x>0.

 SEE ALSO:
 */
{
	// k2 is the value for
	// gamma_R(5./6)*2^(-1./6)
	const double k2 = 1.00563491799858928388289314170833;
	const double k3 = 1.25331413731550012081;   //  sqrt(pi/2)
	const double a1 = 0.22222222222222222222;   //  2/9
	const double a2 = -0.08641975308641974829;  //  -7/89
	const double a3 = 0.08001828989483310284;   // 175/2187
	T res;
	T x_1;

	x_1 = 1. / x;
	res = k2
	      - k3 * exp(-x) * pow(x, (1 / 3.))
    	  * (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)));
	return res;
}
//------------------------------------------------------------------------------------
template<typename T>
inline __device__ double rodconan_gpu(T r, T L0, int k)
/* DOCUMENT rodconan_gpu(r,L0,k=)
 The phase structure function is computed from the expression
 Dphi(r) = k1  * L0^(5./3) * (k2 - (2.pi.r/L0)^5/6 K_{5/6}(2.pi.r/L0))

 For small r, the expression is computed from a development of
 K_5/6 near 0. The value of k2 is not used, as this same value
 appears in the series and cancels with k2.
 For large r, the expression is taken from an asymptotic form.

 SEE ALSO:
 */
{
	const double pi = 3.1415926535897932384626433;
	T res = 0;

	// k1 is the value of :
	// 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
	const double k1 = 0.1716613621245709486;
	const T dprf0 = (2 * pi / L0) * r;
	// k2 is the value for gamma_R(5./6)*2^(-1./6),
	// but is now unused
	// k2 = 1.0056349179985892838;

	// Xlim = 0.75*2*pi;   // = 4.71239
	if (dprf0 > 4.71239)
		res = asymp_macdo_gpu(dprf0);
	else
		res = -macdo_x56_gpu(dprf0, k);

	res *= k1 * pow(L0, 5. / 3);

	return res;
}

template<typename T>
inline __device__ double DPHI_gpu(T x, T y, T indexL0)
/* DOCUMENT dphi = DPHI(x,y,indexL0,rr,tabDPHI,convert) * r0^(-5./3)
 <x> & <y>         :  separation between apertures
 <indexL0>         :  index for the L0 taken into account
 <rr>              :  array of distance between apertures
 NOT USED ANYMORE <tabDPHI>         :  array of precomputed DPHI
 NOT USED ANYMORE <convert>         :  relation between the index on tabDPHI and (x,y)

 Computes the phase structure function for a separation (x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 SEE ALSO:
 */
{
  T r = sqrt(x * x + y * y);

  return rodconan_gpu(r, indexL0, 10);
}


//inline __device__ double cov_XX_gpu(double du, double dv, double ac, double ad, double bc, double bd, double *tabDPHI, double indexL0, double convert, int Ndphi);
//============================================================================================
//============================= MATCOV ELEMENTARY FUNCTIONS ==================================
//============================================================================================
template<typename T>
inline __device__ double cov_XX_gpu(T du, T dv, T ac, T ad, T bc, T bd, T indexL0)
 /* DOCUMENT
   Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return -DPHI_gpu(du + ac, dv, indexL0)
    + DPHI_gpu(du + ad, dv, indexL0)
    + DPHI_gpu(du + bc, dv, indexL0)
    - DPHI_gpu(du + bd, dv, indexL0);;
}

//------------------------------------------------------------------------------------
//inline __device__ double cov_YY_gpu(double du, double dv, double ac, double ad, double bc, double bd, double *tabDPHI, double indexL0, double convert, int Ndphi);
template<typename T>
inline __device__ double cov_YY_gpu(T du, T dv, T ac, T ad, T bc, T bd, T indexL0)
/* DOCUMENT
   Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return  -DPHI_gpu(du, dv + ac, indexL0)
    + DPHI_gpu(du, dv + ad, indexL0)
    + DPHI_gpu(du, dv + bc, indexL0)
    - DPHI_gpu(du, dv + bd, indexL0);;
}

//------------------------------------------------------------------------------------
template<typename T>
inline __device__ double cov_XY_gpu(T du, T dv, T s1, T s2, T indexL0)
{
  return -DPHI_gpu(du + s1 , dv - s2, indexL0)
        + DPHI_gpu(du + s1, dv + s2, indexL0)
        + DPHI_gpu(du - s1, dv - s2, indexL0)
        - DPHI_gpu(du - s1, dv + s2, indexL0);
}

//------------------------------------------------------------------------------------
template<typename T>
inline __device__ double cov_YX_gpu(T du, T dv, T s1, T s2, T indexL0)
{
  return -DPHI_gpu(du - s2, dv + s1, indexL0)
        + DPHI_gpu(du + s2, dv + s1, indexL0)
        + DPHI_gpu(du - s2, dv - s1, indexL0)
        - DPHI_gpu(du + s2, dv - s1, indexL0);
}



//============================================================================================
//=========================== MATCOV 4 (NOISE) KERNELS/FUNCTION ==============================
//============================================================================================
template<typename T>
inline __device__ T compute_element_4(int ipos, int jpos, T *sspSizeL, long *Nssp, T *u, T *v,
					T *indexL0, T *cn2, int Nw, int Nlayer,
					long * Nsubap_wfs, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx,
                                        T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T teldiam)
{
	/* *** Covariance matrix per-element generation ***
	*   Arguments
	*   =========
	*	ipos:		Integer: global x-coordinate of the element w.r.t. the entire matrix
	*	jpos:		Integer: global y-coordinate of the element w.r.t. the entire matrix
	*/

	// for now return a dummy value

  const T lambda2 = 0.00026942094446267851;
  //WFS m
  long Nsubapx = Nsubap_wfs[0];
  int m = 0;
  if (type_mat == 3 || type_mat==4){
       m = Nw-1;
       Nsubapx=Nx;
       ipos+=2*(Nx-Nsubap_wfs[m]);
  }
  else{
      while((ipos / (2 * Nsubapx)) >= 1){
                  m++;
                  Nsubapx += Nsubap_wfs[m];
    }
  }
  Nsubapx -= Nsubap_wfs[m];
  
  //WFS n
  long Nsubapy = Nsubap_wfs[0];
  int n = 0;
  if (type_mat == 2 || type_mat==4){
      n = Nw-1;
       Nsubapy=Nx;
       jpos+=2*(Nx-Nsubap_wfs[m]);
  }
  else{
      while((jpos / (2 * Nsubapy)) >= 1){
        n++;
        Nsubapy += Nsubap_wfs[n];
      }
  }
  Nsubapy -= Nsubap_wfs[n];

  //subap i
  int i = ipos - 2 * Nsubapx;
  //subap j
  int j = jpos - 2 * Nsubapy;
  //xy i
  int xy_i;
  //xy j
  int xy_j;
  if (i>=Nsubap_wfs[m]) {
    i-= Nsubap_wfs[m];
    xy_i = 1;
  } else xy_i = 0;
  if (j>=Nsubap_wfs[n]) {
    j-= Nsubap_wfs[n];
    xy_j = 1;
  } else xy_j = 0;

  const T sspSizem = teldiam / Nssp[m];
  const T sspSizen = teldiam / Nssp[n];

  const T kk = lambda2 / (sspSizem * sspSizen);

  int type = xy_i * 2 + xy_j;

  //Layer l
  T covar = 0.0;
  #pragma unroll
  for (int l = 0; l < Nlayer; l++)
  {
    T sspSizeml = sspSizeL[m * Nlayer + l];
    T sspSizenl = sspSizeL[n * Nlayer + l];
    //test if the altitude layers is not higher than the LGS altitude
    if ((sspSizeml > 0) && (sspSizenl > 0))
    {
      int pos1 = i + Nsubapx + l * Nx;
      int pos2 = j + Nsubapy + l * Nx;
      T du =  u[pos1] - u[pos2];
      T dv =  v[pos1] - v[pos2];

      T s1 = sspSizeml * 0.5;
      T s2 = sspSizenl * 0.5;

      T ac = s1 - s2;
      T ad = s1 + s2;
      T bc = -ad;   // initially -s1-s2;
      T bd = -ac;   // initially -s1+s2;

      if (type == 0){
          covar += 0.5 * cov_XX_gpu(du,dv,ac,ad,bc,bd,indexL0[l]) * kk * cn2[l];
      }
      else if (type == 3){
          covar += 0.5 * cov_YY_gpu(du,dv,ac,ad,bc,bd,indexL0[l]) * kk * cn2[l];
      }
      else if(type==1){
          covar += 0.5 * cov_XY_gpu(du,dv,s1,s2,indexL0[l]) * kk * cn2[l];
      }
      else if(type==2){
          covar += 0.5 * cov_YX_gpu(du,dv,s1,s2,indexL0[l]) * kk * cn2[l];
      }
    }
  }
  // adding noise

  if (m == n && type_mat<4) {
    if (m < nlgs) {
      if (i == j) {
	    // lgs case
	    const int pos1 = i + Nsubapx;
        if (type == 0) covar+=noiseLGSxx[pos1];
        else if (type == 3) covar += noiseLGSyy[pos1];
        else covar += noiseLGSxy[pos1];
      }
      if ((type == 0) || (type == 3))
	covar += lgs_cst;
    } else {
    // ngs case
      if (i==j) {
	if ((type == 0) || (type == 3)) {
	  covar += noiseNGS[m-nlgs];
	}
      }
    }
  }

  return covar;
}





//------------------------------------------------------------------------------------------
template<typename T>
__global__ void matcov_kernel_gpu_4(char uplo, char copy, T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
T *sspSizeL, long *Nssp, T *u, T *v, 
T *L0diff, T *cn2, int Nw, int Nlayer, 
long *Nsubap,long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx, T *noiseLGSyy, T *noiseLGSxy,
int type_mat, int nlgs, T teldiam)
{
  /* *** covariance matrix generation kernel ***
  *  The kernel generates the element values in a given matrix/submatrix
  *   The generation function can be any function, as long as each element 
  *   can be computed both individually and independently
  *
  *  see argument description in the kernel driver
  */
  
  // local thread coordinates w.r.t. thread block
  const int tx_ = threadIdx.x;
  const int ty_ = threadIdx.y;
  
  // local thread block coordinates w.r.t. kernel grid
  const int bx_ = blockIdx.x;
  const int by_ = blockIdx.y;
  
  // local coordinates of the element w.r.t. submatrix
  int lx = bx_ * blockDim.x + tx_;
  int ly = by_ * blockDim.y + ty_;
  
  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;
  
  // out-of-bound threads should terminate
  if( (lx >= nrows) || (ly >= ncols) ) return;
  
  // Advance the data pointer accordingly
  //data += ly * lda + lx;
  
  T value;
  if(uplo == 'l')
  {
      if(gy <= gx)
      {   
          value = compute_element_4(gx, gy, sspSizeL, Nssp, u, v,
                                    L0diff, cn2, Nw, Nlayer, Nsubap, 
                                    Nx, lgs_cst, noiseNGS,noiseLGSxx,noiseLGSyy,noiseLGSxy,
                                    type_mat,nlgs,teldiam);
          data[ly * lda + lx] = value;
          if(copy == 'c') data[lx * lda + ly] = value;
      }
  }
  else if (uplo == 'u') // upper
  {
      if(gx <= gy)
      {
          value = compute_element_4(gx, gy, sspSizeL, Nssp, u, v,
                                    L0diff, cn2, Nw, Nlayer, Nsubap, 
                                    Nx, lgs_cst, noiseNGS,noiseLGSxx,noiseLGSyy,noiseLGSxy ,
                                    type_mat,nlgs,teldiam);
          data[ly * lda + lx] = value;
          if(copy == 'c') data[lx * lda + ly] = value;
      }
  }
  else    // uplo = 'f' full generation
  {
          value = compute_element_4(gx, gy, sspSizeL, Nssp, u, v,
                                    L0diff, cn2, Nw, Nlayer, Nsubap, 
                                    Nx, lgs_cst, noiseNGS,noiseLGSxx,noiseLGSyy,noiseLGSxy ,
                                    type_mat,nlgs,teldiam);
      data[ly * lda + lx] = value;
  }
  
}

template __global__ void matcov_kernel_gpu_4(char uplo, char copy, float* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
float *sspSizeL, long *Nssp, float *u, float *v, 
float *L0diff, float *cn2, int Nw, int Nlayer, 
long *Nsubap,long Nx, float lgs_cst, float *noiseNGS, float *noiseLGSxx, float *noiseLGSyy, float *noiseLGSxy,
int type_mat, int nlgs, float teldiam);
template __global__ void matcov_kernel_gpu_4(char uplo, char copy, double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
double *sspSizeL, long *Nssp, double *u, double *v, 
double *L0diff, double *cn2, int Nw, int Nlayer, 
long *Nsubap,long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx, double *noiseLGSyy, double *noiseLGSxy,
int type_mat, int nlgs, double teldiam);









//============================================================================================
//============================= MATCOV 4 (NOISE) =============================================
//============================================================================================
#include <curand.h>
template<typename T>
void matcov_tile_gpu_4(  char uplo, char copy, T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  T *sspSizeL, long *Nssp, T *u, T *v,
  T *indexL0, T *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx,
  T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T DiamTel,cudaStream_t stream)
{
	/* *** matcov gpu kernel driver ***
	*  Arguments
	*  ==========
	*  data		double pointer: A pointer to the matrix/submatrix to be generated. It  
	*  			should always point to the first element in a matrix/submatrix
	*
	*  nrows	integer: The number of rows of the matrix/submatrix to be generated	
	*
	*  ncols	integer: The number of columns of the matrix/submatrix to be generated
	*
	*  xoffset	integer: The x-offset of the submatrix, must be zero if the entire matrix
	*			is generated. Its the x-coordinate of the first element in the matrix/submatrix
	*
	*  yoffset  integer: The y-offset of the submatrix, must be zero if the entire matrix
	*			is generated. Its the y-coordinate of the first element in the matrix/submatrix
	*
	*  lda		integer: The leading dimension of the matrix/submatrix
	*/

  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  

  matcov_kernel_gpu_4<T><<<dimGrid, dimBlock, 0, stream>>>(uplo, copy, data, nrows, ncols, xoffset, yoffset, lda, sspSizeL, 
						Nssp, u, v,
						indexL0, cn2, Nw, Nlayer, 
						Nsubap, Nx, lgs_cst, noiseNGS, noiseLGSxx,
						noiseLGSyy, noiseLGSxy,type_mat, nlgs, DiamTel);

	cudaStreamSynchronize(stream);
}

template void matcov_tile_gpu_4<float>(  char uplo, char copy, float* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  float *sspSizeL, long *Nssp, float *u, float *v,
  float *indexL0, float *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, float lgs_cst, float *noiseNGS, float *noiseLGSxx,
  float *noiseLGSyy, float *noiseLGSxy, int type_mat, int nlgs, float DiamTel,cudaStream_t stream);

template void matcov_tile_gpu_4<double>(  char uplo, char copy, double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  double *sspSizeL, long *Nssp, double *u, double *v,
  double *indexL0, double *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx,
  double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double DiamTel,cudaStream_t stream);
