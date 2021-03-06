#include<cuda_runtime.h>
#include<cuda_runtime_api.h>
#include<cublas.h>
#include<stdio.h>
#include "matcov.h"
#include "matcov_gpu_gb.h"

/*  Tuning parameters of tbulateDPHI kernel*/
#define tabDPHI_thread_x	(256)

/*	Tuning parameters of matcov GPU Kernel */
// Thread block size (x, y), 
// max #threads per block is 512 for fermi and 1024 for kepler
#define matcov_thread_x	(8)
#define matcov_thread_y	(8)

//#define CUDA_ERROR_CHECK
 
#define CudaSafeCall( err ) __cudaSafeCall( err, __FILE__, __LINE__ )
#define CudaCheckError()    __cudaCheckError( __FILE__, __LINE__ )
 
inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaSafeCall() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
 
    return;
}
 
inline void __cudaCheckError( const char *file, const int line )
{
#ifdef CUDA_ERROR_CHECK
    cudaError err = cudaGetLastError();
    if ( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
 
    // More careful checking. However, this will affect performance.
    // Comment away if needed.
    err = cudaDeviceSynchronize();
    if( cudaSuccess != err )
    {
        fprintf( stderr, "cudaCheckError() with sync failed at %s:%i : %s\n",
                 file, line, cudaGetErrorString( err ) );
        exit( -1 );
    }
#endif
 
    return;
}

//============================================================================================
//================================= AUX FUNCTIONS ============================================
//============================================================================================
#define VERBOSE 0
void process_err(cudaError_t e, const char* str)
{
	if(VERBOSE) printf("%s\n", str);
	if(e != cudaSuccess)
	{
		printf("*** Error %s: %s \n", str, cudaGetErrorString(e));
		exit(1);
	}
}
//-----------------------------------------------------------------------
double* arr2dAlloc_gpu_gb(long nbLin, long nbCol)
/* DOCUMENT  array = arr2dAlloc(nblin,nbcol)

 Allocates a 2d array (double).
 */
{
	cudaError_t e;
	double* tableau;
	e = cudaMalloc((void**)&tableau, sizeof(double) * nbCol * nbLin);
	process_err(e, "gpu alloc tableau2");
	return tableau;
}

void arr2dFree_gpu_gb(double *tableau)
/* DOCUMENT  arr2dFree(array)

 Free a 2d array (double).
 */
{
	if(tableau)cudaFree(tableau);
}

//============================================================================================
//============================= tabDPHI FUNCTIONS/KERNEL(s) ==================================
//============================================================================================
__device__ double macdo_x56_gpu_gb(double x, int k)
/* DOCUMENT  macdo_x56_gpu_gb(x)

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
	double x2n;               // x^2.a, etc
	double s = 0.0;
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
__device__ double asymp_macdo_gpu_gb(double x)
/* DOCUMENT asymp_macdo_gpu_gb(x)

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
	double res;
	double x_1;
	
	x_1 = 1. / x;
	res = k2
	      - k3 * exp(-x) * pow(x, 1 / 3.)
    	  * (1.0 + x_1 * (a1 + x_1 * (a2 + x_1 * a3)));
	return res;
}
//------------------------------------------------------------------------------------
__device__ double rodconan_gpu_gb(double r, double L0, int k)
/* DOCUMENT rodconan_gpu_gb(r,L0,k=)
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
	double res = 0;
	
	// k1 is the value of :
	// 2*gamma_R(11./6)*2^(-5./6)*pi^(-8./3)*(24*gamma_R(6./5)/5.)^(5./6);
	const double k1 = 0.1716613621245709486;
	const double dprf0 = (2 * pi / L0) * r;
	// k2 is the value for gamma_R(5./6)*2^(-1./6),
	// but is now unused
	// k2 = 1.0056349179985892838;
	
	// Xlim = 0.75*2*pi;   // = 4.71239
	if (dprf0 > 4.71239)
		res = asymp_macdo_gpu_gb(dprf0);
	else
		res = -macdo_x56_gpu_gb(dprf0, k);

	res *= k1 * pow(L0, 5. / 3);
	return res;
}

__global__ void tabulateDPHI_gpu_gb_kernel(double* tabDPHI_d, double* L0diff_d, long Nl0, long Ndphi, double convert)
{
	const int tx = threadIdx.x;
	const int ty = blockIdx.x;
	
	const int tid = ty * blockDim.x + tx;
	int l = tid / Ndphi;
	int j = tid % Ndphi;
	
	if(tid >= (Nl0*Ndphi) ) return;
	
	tabDPHI_d[tid] = rodconan_gpu_gb((double)j / convert, L0diff_d[l], 10);
	
	//double* mytabDPHI = tabDPHI_d + (l * Ndphi);
	//
	//int j, k;
	//#pragma unroll
	//for(k = 0; k < (Ndphi/tabDPHI_thread_x); k++)
	//{
	//	j = k * tabDPHI_thread_x + tx;
	//	mytabDPHI[j] = rodconan_gpu_gb(rr_d[j], L0diff_d[l], 10);
	//}
	//
	//k = (Ndphi/tabDPHI_thread_x);
	//if(tx < (Ndphi%tabDPHI_thread_x) )
	//{
	//	j = k * tabDPHI_thread_x + tx;
	//	mytabDPHI[j] = rodconan_gpu_gb(rr_d[j], L0diff_d[l], 10);
	//}
}
//------------------------------------------------------------------------------------
__device__ double DPHI_gpu_gb(double x, double y, long indexL0, double *tabDPHI, double convert, int Ndphi)
/* DOCUMENT dphi = DPHI(x,y,indexL0,rr,tabDPHI,convert) * r0^(-5./3)
 <x> & <y>         :  separation between apertures
 <indexL0>         :  index for the L0 taken into account
 <rr>              :  array of distance between apertures
 <tabDPHI>         :  array of precomputed DPHI
 <convert>         :  relation between the index on tabDPHI and (x,y)

 Computes the phase structure function for a separation (x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 SEE ALSO:
 */
{
  double r = sqrt(x * x + y * y);
  long i0 = (long) (r * convert);
  long i1 = i0 + 1;

  return ((r - (double)i0 / convert) * tabDPHI[indexL0 * Ndphi + i1]
	  + ((double)i1 / convert - r) * tabDPHI[indexL0 * Ndphi + i0]);
      
}

//============================================================================================
//============================= SUBAP POSITION KERNELS/FUNCTIONS =============================
//============================================================================================
__global__ void subposition_gpu_gb_kernel(long Nw, long Nsubap, long Nlayer, double *alphaX, double *alphaY,
				       double *h, double *GsAlt, long *Nssp, double *diamPup, double *thetaML,
				       long *ioff, double *X, double *Y, double *XPup, double *YPup, 
				       double *u, double *v)
{
  const int tx = threadIdx.x;
  const int ty = blockIdx.x;
	
  const int tid = ty * blockDim.x + tx;
  long i;
  long n;
  long l;
  const double rad = 3.14159265358979323846 / 180.;

  if(tid >= (Nw * Nsubap * Nlayer) ) return;
	
  l = tid / (Nw * Nsubap);

  const int pos = tid - l * (Nsubap * Nw);

  i = pos / Nw;
  n = pos - i * Nw;

  //tid = n + i * Nw + l * Nw * Nsubap

  const double dX = alphaX[n] * h[l];
  const double dY = alphaY[n] * h[l];
    
  const double rr = 1. - h[l] * GsAlt[n];
    
  const long nssp = Nssp[n];
    
  //magnification factor
  const double G = diamPup[n] / (double) (nssp);
    
  //rotation angle
  const double th = thetaML[n] * rad;

  //taking magnification factor into account
  const double xtp = X[ioff[n] + i] * G;
  const double ytp = Y[ioff[n] + i] * G;
    
  //taking rotation into account
  double uu = xtp * cos(th) - ytp * sin(th);
  double vv = xtp * sin(th) + ytp * cos(th);
    
  //taking pupil offset into account
  uu += XPup[n];
  vv += YPup[n];
    
  //Projection onto  the layer
  u[tid] = uu * rr + dX;
  v[tid] = vv * rr + dY;
}

//============================================================================================
//============================= MATCOV ELEMENTARY FUNCTIONS ==================================
//============================================================================================
__device__ double cov_XX_gpu_gb(double du, double dv, double ac, double ad, double bc, double bd, double *tabDPHI, long indexL0, double convert, int Ndphi)
 /* DOCUMENT
   Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return -DPHI_gpu_gb(du + ac, dv, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du + ad, dv, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du + bc, dv, indexL0, tabDPHI, convert, Ndphi)
    - DPHI_gpu_gb(du + bd, dv, indexL0, tabDPHI, convert, Ndphi);
}

//------------------------------------------------------------------------------------
__device__ double cov_YY_gpu_gb(double du, double dv, double ac, double ad, double bc, double bd, double *tabDPHI, long indexL0, double convert, int Ndphi)
/* DOCUMENT
   Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{ 
  return  -DPHI_gpu_gb(du, dv + ac, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du, dv + ad, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du, dv + bc, indexL0, tabDPHI, convert, Ndphi)
    - DPHI_gpu_gb(du, dv + bd, indexL0, tabDPHI, convert, Ndphi);
}


//------------------------------------------------------------------------------------
__device__ double cov_XY_gpu_gb(double du, double dv, double s0, double *tabDPHI, long indexL0, double convert, int Ndphi)
/* DOCUMENT
   Compute the XY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
{
  return -DPHI_gpu_gb(du + s0, dv - s0, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du + s0, dv + s0, indexL0, tabDPHI, convert, Ndphi)
    + DPHI_gpu_gb(du - s0, dv - s0, indexL0, tabDPHI, convert, Ndphi)
    - DPHI_gpu_gb(du - s0, dv + s0, indexL0, tabDPHI, convert, Ndphi);
}


//============================================================================================
//============================= MATCOV 3 FUNCTIONS/KERNEL ====================================
//============================================================================================
__device__ double compute_element_3(int ipos, int jpos, double convert,
				  double *sspSizeL, long *Nssp, double *u, double *v, double pasDPHI,double *tabDPHI, 
				  long *indexL0, double *cn2, int Ndphi, int Nw, int Nlayer, int Nsubap,
				  int type_mat, double teldiam)
{
	/* *** Covariance matrix per-element generation ***
	*   Arguments
	*   =========
	*	ipos:		Integer: global x-coordinate of the element w.r.t. the entire matrix
	*	jpos:		Integer: global y-coordinate of the element w.r.t. the entire matrix
	*/
	
  const double lambda2 = 0.00026942094446267851;
  const int nslps = Nsubap*2;
  
  //WFS m
  int m = ipos / nslps; //tab_wfs[ipos];
  if (type_mat == 3) m = Nw-1;
  //WFS n
  int n = jpos / nslps; //tab_wfs[jpos];
  if (type_mat == 2) n = Nw-1;
  
  //subap i
  int i = ipos % (nslps/2); //tab_subap[ipos];
  //subap j
  int j = jpos % (nslps/2); //tab_subap[jpos];;
  
  //xy i
  int xy_i = (ipos / (nslps/2))%2;  //tab_xy[ipos]; 
  //xy j
  int xy_j = (jpos / (nslps/2))%2;  //tab_xy[jpos];
  
  const double sspSizem = teldiam / Nssp[m];
  const double sspSizen = teldiam / Nssp[n];
  
  const double kk = lambda2 / (sspSizem * sspSizen);
    
  int type = xy_i * 2 + xy_j;

  //Layer l
  double covar = 0.0;
  #pragma unroll
  for (int l = 0; l < Nlayer; l++) 
  {
    const double sspSizeml = sspSizeL[m * Nlayer + l];
    const double sspSizenl = sspSizeL[n * Nlayer + l];
    //test if the altitude layers is not higher than the LGS altitude
    if ((sspSizeml > 0) && (sspSizenl > 0)) 
    {
      const int pos1 = m + i * Nw + l * Nw * Nsubap;
      const int pos2 = n + j * Nw + l * Nw * Nsubap;
      const double du = u[pos1] - u[pos2];	      
      const double dv =  v[pos1] - v[pos2];
      
      const double s1 = sspSizeml * 0.5;
      const double s2 = sspSizenl * 0.5;
      
      const double ac = s1 - s2;
      const double ad = s1 + s2;
      const double bc = -ad;   // initially -s1-s2;
      const double bd = -ac;   // initially -s1+s2;

      if (type == 0) covar += 0.5 * pasDPHI * cov_XX_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else if (type == 3) covar += 0.5 * pasDPHI * cov_YY_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else //if ((type == 1) || (type == 2)) 
      {
      	const double s0 = sqrt(s1 * s1 + s2 * s2); //half size of the subaperture equivalent to a convolution by s1 and s2
      	const double dd = (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2; // Nono's style ....
      	covar += 0.25 * pasDPHI * cov_XY_gpu_gb(du,dv,s0,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  return (double)covar; 
}

__global__ void matcov_kernel_3(char uplo, char copy, double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
					double convert, double *sspSizeL, long *Nssp, double *u, double *v, double pasDPHI,double *tabDPHI, long *indexL0, 
				  	double *cn2, int Ndphi, int Nw, int Nlayer, int Nsubap, int type_mat, double teldiam)
{
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element 
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
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
	
  double value;
  if(uplo == 'l')
  {
  	if(gy <= gx)
  	{	
  		value = compute_element_3(gx, gy,convert,sspSizeL,Nssp,u,v,pasDPHI,tabDPHI, indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,type_mat,teldiam);
  		data[ly * lda + lx] = value;
  		if(copy == 'c') data[lx * lda + ly] = value;
  	}
  }
  else if (uplo == 'u') // upper
  {
  	if(gx <= gy)
  	{
  		value = compute_element_3(gx, gy,convert,sspSizeL,Nssp,u,v,pasDPHI,tabDPHI, indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,type_mat,teldiam);
  		data[ly * lda + lx] = value;
  		if(copy == 'c') data[lx * lda + ly] = value;
  	}
  }
  else	// uplo = 'f' full generation
  {
  	value = compute_element_3(gx, gy,convert,sspSizeL,Nssp,u,v,pasDPHI,tabDPHI, indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,type_mat,teldiam);
  	data[ly * lda + lx] = value;
  }
  	
  //if ((type_mat == 3) || (gx <= gy)) 
  //{
    // call the generation function
    //data[0] = compute_element_3(gx, gy, tab_wfs, tab_subap, tab_xy,convert,sspSizeL,Nssp,u,v,pasDPHI,tabDPHI,
	//		      indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,type_mat,teldiam);
    //printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
  //} 
}

//============================================================================================
//============================= MATCOV TS FUNCTIONS/KERNEL ===================================
//============================================================================================
__device__ double compute_element_ts_(int ipos, int jpos, double convert, double *X, double *Y, 
				     long *Nssp, double pasDPHI, double *tabDPHI, long *indexL0, double *cn2, 
				     int Ndphi, int Nw, int Nlayer, int Nsubap, double teldiam)
{
	/* *** Covariance matrix per-element generation ***
	*   Arguments
	*   =========
	*	ipos:		Integer: global x-coordinate of the element w.r.t. the entire matrix
	*	jpos:		Integer: global y-coordinate of the element w.r.t. the entire matrix
	*/
	
	// for now return a dummy value
  
  double lambda2 = 0.00026942094446267851;
  //WFS Nw-1
   //subap i
  int i = ipos < Nsubap ? ipos : ipos - Nsubap;
  //subap j
  int j = jpos < Nsubap ? jpos : jpos - Nsubap;
  //xy i
  int xy_i = ipos < Nsubap ? 0 : 1;
  //xy j
  int xy_j = jpos < Nsubap ? 0 : 1;
  
  double sspSize = teldiam / Nssp[Nw-1];
  
  double kk = lambda2 / (sspSize * sspSize);
    
  int type = xy_i * 2 + xy_j;

  double s = sspSize * 0.5;
  
  double ac = 0.0;
  double ad = 2.0 * s;
  double bc = -ad;   
  double bd = 0.0;   

  double du = X[(Nsubap*(Nw-1)+i)] - X[(Nsubap*(Nw-1)+j)];	      
  double dv = Y[(Nsubap*(Nw-1)+i)] - Y[(Nsubap*(Nw-1)+j)];
  
  //if(ipos < 10)printf("ipos = %d - %d\n", ipos, (Nsubap*(Nw-1)+i));
  //if(jpos < 10)printf("jpos = %d - %d\n", jpos, (Nsubap*(Nw-1)+j));
  
  //const double du = X[0] - X[1];	      
  //const double dv = Y[0] - Y[1];
  
//Layer l
  double covar = 0.0;
  #pragma unroll
  for (int l = 0; l < Nlayer; l++) 
  {
     //test if the altitude layers is not higher than the LGS altitude
    if (sspSize > 0) 
    {
      if (type == 0) covar += 0.5 * pasDPHI * cov_XX_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else if (type == 3) covar += 0.5 * pasDPHI * cov_YY_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else 
      {
      	double s0 = 1.41421*s; //half size of the subaperture equivalent to a convolution by s1 and s2
      	double dd = 0;
      	covar += 0.25 * pasDPHI * cov_XY_gpu_gb(du,dv,s0,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  return (double)covar; 
}
//--------------------------------------------------------------------------------------------
__global__ void matcov_ts_kernel(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
				  double convert, double *X, double *Y, long *Nssp, double pasDPHI,double *tabDPHI, 
				 long *indexL0, double *cn2, int Ndphi, int Nw, int Nlayer, int Nsubap, double teldiam)
{
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element 
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
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
  data += ly * lda + lx;
	
    // call the generation function
    data[0] = compute_element_ts_(gx, gy, convert,X, Y,Nssp,pasDPHI,tabDPHI, indexL0,cn2,Ndphi,Nw,Nlayer,Nsubap,teldiam);
    //printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
}

//============================================================================================
//============================= MATCOV TS  ===================================================
//============================================================================================
void matts_gpu_gb(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu)
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
	
  const long Nw = tomo.Nw;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);


  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = tomo.Nsubap[Nw-1];
  
  matcov_ts_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(data, nrows, ncols, xoffset, yoffset, lda, 
					   convert,tomo_gpu->X_d,tomo_gpu->Y_d,tomo_gpu->Nssp_d,
					   pasDPHI,tomo_gpu->tabDPHI_d,tomo_gpu->indexL0_d,tomo_gpu->cn2_d,
					   Ndphi,tomo.Nw,tomo.Nlayer,Nsubap,tomo.DiamTel);
  //CudaCheckError();
}

//============================================================================================
//============================= MATCOV COPY KERNEL ===========================================
//============================================================================================
__global__ void matcov_kernel_copy(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda)
{
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element 
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
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
	
  if (gx > gy) {
    // call the generation function
    data[ly * lda + lx] = data[ly + lx * lda];
    //printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
  }
}

//============================================================================================
//============================= MATCOV 1 =====================================================
//============================================================================================
//************************** OBSOLETE - REMOVED ********************************************//

//============================================================================================
//============================= MATCOV 2 =====================================================
//============================================================================================
//************************** OBSOLETE - REMOVED ********************************************//

//============================================================================================
//=============================== TOMO INIT/FIN FUNCTIONS ====================================
//============================================================================================
void init_tomo_gpu_gb(struct gtomo_struct *tomo_gpu, struct tomo_struct tomo)
{
  cudaError_t e;

  e = cudaMalloc((void**)&(tomo_gpu->indexL0_d), tomo.Nlayer*sizeof(long));
  process_err(e, "alloc gpu indexL0_d");

  e = cudaMalloc((void**)&(tomo_gpu->u_d), tomo.Nlayer*tomo.Nsubap[0]*tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu u_d");
  //printf("size of u is %d\n", tomo.Nlayer*tomo.Nsubap[0]*tomo.Nw);
  //printf("u_d = 0x%x \n", (tomo_gpu->u_d) );
  
  e = cudaMalloc((void**)&(tomo_gpu->v_d), tomo.Nlayer*tomo.Nsubap[0]*tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu v_d");
  //printf("size of v is %d\n", tomo.Nlayer*tomo.Nsubap[0]*tomo.Nw);
  //printf("v_d = 0x%x \n", (tomo_gpu->v_d) );
  
  e = cudaMalloc((void**)&(tomo_gpu->sspSizeL_d), tomo.Nw*tomo.Nlayer*sizeof(double));
  process_err(e, "alloc gpu sspSizeL_d");

  e = cudaMalloc((void**)&(tomo_gpu->cn2_d), tomo.Nw*tomo.Nlayer*sizeof(double));
  process_err(e, "alloc gpu cn2_d");

  e = cudaMalloc((void**)&(tomo_gpu->h_d), tomo.Nlayer*sizeof(double));
  process_err(e, "alloc gpu h_d");


  e = cudaMalloc((void**)&(tomo_gpu->Nssp_d), tomo.Nw*sizeof(long));
  process_err(e, "alloc gpu Nssp_d");

  e = cudaMalloc((void**)&(tomo_gpu->ioff_d), tomo.Nw*sizeof(long));
  process_err(e, "alloc gpu ioff_d");

  e = cudaMalloc((void**)&(tomo_gpu->alphaX_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu alphaX_d");
  
  e = cudaMalloc((void**)&(tomo_gpu->alphaY_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu alphaY_d");
  
  e = cudaMalloc((void**)&(tomo_gpu->GsAlt_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu GsAlt_d");

  e = cudaMalloc((void**)&(tomo_gpu->diamPup_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu diamPup_d");

  e = cudaMalloc((void**)&(tomo_gpu->thetaML_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu thetaML_d");

  e = cudaMalloc((void**)&(tomo_gpu->X_d), tomo.Nx*sizeof(double));
  process_err(e, "alloc gpu X_d");
  //printf("size of X is %d\n", tomo.Nx);
  //printf("X_d = 0x%x \n", (tomo_gpu->X_d) );
  
  e = cudaMalloc((void**)&(tomo_gpu->Y_d), tomo.Nx*sizeof(double));
  process_err(e, "alloc gpu Y_d");
  //printf("size of X is %d\n", tomo.Nx);
  //printf("Y_d = 0x%x \n", (tomo_gpu->Y_d) );
  
  e = cudaMalloc((void**)&(tomo_gpu->XPup_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu XPup_d");

  e = cudaMalloc((void**)&(tomo_gpu->YPup_d), tomo.Nw*sizeof(double));
  process_err(e, "alloc gpu YPup_d");
  /*
  e = cudaMalloc((void**)&(tomo_gpu->Cmm_d), tomo.Nw*tomo.Nsubap[0]*2*tomo.Nw*tomo.Nsubap[0]*2*sizeof(double));
  process_err(e, "alloc gpu YPup_d");

  e = cudaMalloc((void**)&(tomo_gpu->Cpm_d), tomo.Nsubap[0]*2*tomo.Nw*tomo.Nsubap[0]*2*sizeof(double));
  process_err(e, "alloc gpu YPup_d");

  e = cudaMalloc((void**)&(tomo_gpu->R_d), tomo.Nsubap[0]*2*tomo.Nw*tomo.Nsubap[0]*2*sizeof(double));
  process_err(e, "alloc gpu YPup_d");
  */
  
  tomo_gpu->L0diff_d = NULL;
  tomo_gpu->tabDPHI_d = NULL;
  
  e = cudaStreamCreate(&(tomo_gpu->matcov_stream));
  process_err(e, "create matcov stream");

}

void free_tomo_gpu_gb(struct gtomo_struct *tomo_gpu)
{
  cudaError_t e;

  if ((tomo_gpu->u_d) != NULL) e = cudaFree(tomo_gpu->u_d);
  process_err(e, "free gpu u_d");

  if (tomo_gpu->v_d) e = cudaFree(tomo_gpu->v_d);
  process_err(e, "free gpu v_d");

  if (tomo_gpu->sspSizeL_d) e = cudaFree(tomo_gpu->sspSizeL_d) ;
  process_err(e, "free gpu sspSizeL_d");

  if (tomo_gpu->cn2_d) e = cudaFree(tomo_gpu->cn2_d);
  process_err(e, "free gpu cn2_d");

  if (tomo_gpu->h_d) e = cudaFree(tomo_gpu->h_d);
  process_err(e, "free gpu h_d");

  if (tomo_gpu->indexL0_d) e = cudaFree(tomo_gpu->indexL0_d);
  process_err(e, "free gpu indexL0_d");


  if (tomo_gpu->Nssp_d) e = cudaFree(tomo_gpu->Nssp_d);
  process_err(e, "free gpu Nssp_d");

  if (tomo_gpu->ioff_d) e = cudaFree(tomo_gpu->ioff_d);
  process_err(e, "free gpu ioff_d");

  if (tomo_gpu->alphaX_d) e = cudaFree(tomo_gpu->alphaX_d);
  process_err(e, "free gpu alphaX_d");

  if (tomo_gpu->alphaY_d) e = cudaFree(tomo_gpu->alphaY_d);
  process_err(e, "free gpu alphaY_d");

  if (tomo_gpu->GsAlt_d) e = cudaFree(tomo_gpu->GsAlt_d);
  process_err(e, "free gpu GsAlt_d");

  if (tomo_gpu->diamPup_d) e = cudaFree(tomo_gpu->diamPup_d);
  process_err(e, "free gpu diamPup_d");

  if (tomo_gpu->thetaML_d) e = cudaFree(tomo_gpu->thetaML_d);
  process_err(e, "free gpu thetaML_d");

  if (tomo_gpu->X_d) e = cudaFree(tomo_gpu->X_d);
  process_err(e, "free gpu X_d");

  if (tomo_gpu->Y_d) e = cudaFree(tomo_gpu->Y_d);
  process_err(e, "free gpu Y_d");

  if (tomo_gpu->XPup_d) e = cudaFree(tomo_gpu->XPup_d);
  process_err(e, "free gpu XPup_d");

  if (tomo_gpu->YPup_d) e = cudaFree(tomo_gpu->YPup_d);
  process_err(e, "free gpu YPup_d");

  /*
  if (tomo_gpu->Cmm_d) e = cudaFree(tomo_gpu->Cmm_d);
  process_err(e, "free gpu YPup_d");

  if (tomo_gpu->Cpm_d) e = cudaFree(tomo_gpu->Cpm_d);
  process_err(e, "free gpu YPup_d");

  if (tomo_gpu->R_d) e = cudaFree(tomo_gpu->R_d);
  process_err(e, "free gpu YPup_d");
  */

  if ((tomo_gpu->tabDPHI_d) != NULL) e = cudaFree(tomo_gpu->tabDPHI_d);
  process_err(e, "free gpu tabDPHI_d");

  if ((tomo_gpu->L0diff_d) != NULL) e = cudaFree(tomo_gpu->L0diff_d);
  process_err(e, "free gpu L0diff_d");
  
  // destroy matcov stream
  e = cudaStreamDestroy(tomo_gpu->matcov_stream);
  process_err(e, "destroy matcov stream");
}


//============================================================================================
//============================ MATCOV V3/V4 DPHI/SUBAP FUNCTIONS =============================
//============================================================================================
void tab_dphi_gpu_gb(double *tab_dphi, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu, long Ndphi, double *L0diff_d, int Nl0, double convert)
//void tabulateDPHI_gpu_gb(double* tabDPHI_d, double* rr_d,struct tomo_struct tomo, long Ndphi, long *indexL0_h)
/* DOCUMENT tabDPHI = tabulateDPHI(rr,tomo,Ndphi, indexL0)
 <tomo>            :  structure with all the needed information
 <Ndphi>           :  size of rr
 <indexL0>         :  link between the index of the studied layer and the index of the precomputed one. 

 Computes the phase structure function for a separation rr(x,y).
 The r0 is not taken into account : the final result of DPHI(x,y,L0)
 has to be scaled with r0^-5/3, with r0 expressed in meters, to get
 the right value.

 Computes the phase structure for each different L0 and give a array (indexL0) to link the index of the layer i and the index of tabDPHI : for the layer l, DPHI = DPHI( du, dv, indexL0[l],rr,tabDPHI, convert).
 SEE ALSO: DPHI
 */
{
  // Assume one thread per element
  int nblocks = (Ndphi*Nl0)/tabDPHI_thread_x + ( ((Ndphi*Nl0)%tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);

  tabulateDPHI_gpu_gb_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(tab_dphi, L0diff_d, Nl0, Ndphi, convert);
  //CudaCheckError();
}
//------------------------------------------------------------------------------------
//extern "C"
void sub_pos_gpu_gb(struct gtomo_struct *tomo_gpu, struct tomo_struct tomo)
//void subap_position_gpu_gb(struct tomo_struct tomo, double ***u, double ***v)
/* DOCUMENT DOCUMENT         subap_position(tomo, u, v)
   <tomo>                : structure with all the needed information.
   <u> and <v>           : 3d arrays containing the sub-apertures projected coordinates onto all the layers. u[0][2][1] is the X-coordinate of the subap 2 of the WFS 0 on the layer 1.

   Computes the projected coordinates of all subapertures  projected onto all the layer
 */
{
  int msize = tomo.Nlayer * tomo.Nw * tomo.Nsubap[0];
  int nblocks = msize / tabDPHI_thread_x + ( ( msize % tabDPHI_thread_x) != 0);
  dim3 dimBlock(tabDPHI_thread_x, 1);
  dim3 dimGrid(nblocks, 1);
  subposition_gpu_gb_kernel<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(tomo.Nw, tomo.Nsubap[0], tomo.Nlayer, tomo_gpu->alphaX_d, 
						tomo_gpu->alphaY_d,tomo_gpu->h_d, tomo_gpu->GsAlt_d, 
						tomo_gpu->Nssp_d, tomo_gpu->diamPup_d, tomo_gpu->thetaML_d, 
						tomo_gpu->ioff_d, tomo_gpu->X_d, tomo_gpu->Y_d, 
						tomo_gpu->XPup_d, tomo_gpu->YPup_d, tomo_gpu->u_d, tomo_gpu->v_d);
  //CudaCheckError();

}

//============================================================================================
//=============================== TOMO UPDATE FUNCTIONS ======================================
//============================================================================================
void update_tomo_atm_gpu_gb(struct gtomo_struct *tomo_gpu,struct tomo_struct tomo) 
{

  cudaError_t e;

  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);

  e = cudaMemcpyAsync(tomo_gpu->h_d, tomo.h, tomo.Nlayer*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu h_d");
  e = cudaMemcpyAsync(tomo_gpu->cn2_d, tomo.cn2, tomo.Nlayer*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu cn2_d");


  double *sspSizeL = (double *)malloc(sizeof(double)*tomo.Nw*tomo.Nlayer);
  for (int cc = 0; cc < tomo.Nw * tomo.Nlayer; cc++) {
    int n = cc / tomo.Nlayer;
    int l = cc - n * tomo.Nlayer;
    sspSizeL[cc] = tomo.sspSize[n] * (1. - tomo.GsAlt[n] * tomo.h[l]);
  }

  e = cudaMemcpyAsync(tomo_gpu->sspSizeL_d, sspSizeL, tomo.Nw*tomo.Nlayer*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu sspSizeL_d");
  cudaStreamSynchronize(tomo_gpu->matcov_stream);
	
  //Search the different L0 and build indexL0
  const long Nlayer = tomo.Nlayer;
  long i, j;
  int cpt = 1;
  double tmp[Nlayer];
  long indexL0[Nlayer];

  tmp[0] = tomo.L0[0];
  indexL0[0] = 0;
  
  for (i = 1; i < Nlayer; i++) {
    j = 0;
    const double l0 = tomo.L0[i];
    
    while ((j < cpt) && (tmp[j] != l0)) {j++;}
    
    indexL0[i] = j;
    
    if (j == cpt) {
      tmp[j] = l0;
      cpt++;
    }
  }
  
  e = cudaMemcpyAsync((tomo_gpu->indexL0_d), indexL0, tomo.Nlayer*sizeof(long), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu indexL0_d");

  tomo_gpu->Nl0 = cpt;
  double L0diff[tomo_gpu->Nl0];

  // allocate space for L0
  if ((tomo_gpu->L0diff_d) != NULL){cudaFree(tomo_gpu->L0diff_d);}
  e = cudaMalloc((void**)&(tomo_gpu->L0diff_d), tomo_gpu->Nl0*sizeof(double));
  process_err(e, "alloc gpu L0diff_d");
  
  for (i = 0; i < tomo_gpu->Nl0; i++)  {
    L0diff[i] = tmp[i];
  }
  
  // offload L0diff
  e = cudaMemcpyAsync(tomo_gpu->L0diff_d, L0diff, tomo_gpu->Nl0*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "offload L0diff");
  
  //précalcul de DPHI : que pour chaque différent L0
  if ((tomo_gpu->tabDPHI_d) != NULL){cudaFree(tomo_gpu->tabDPHI_d);}
  //printf("tabDPHI alloc \n");
  e = cudaMalloc((void**)&(tomo_gpu->tabDPHI_d), tomo_gpu->Nl0*Ndphi*sizeof(double));
  process_err(e, "alloc gpu tabDPHI_d");
  
  tab_dphi_gpu_gb(tomo_gpu->tabDPHI_d, tomo, tomo_gpu, Ndphi, tomo_gpu->L0diff_d, tomo_gpu->Nl0,convert);

  // %%%%%%% Computation of the sub-apertures positions and sizes %%%%%%%%%%%
 // u, v :arrays containing all the sub-apertures coordinates of all WFS, one after the other
  // u[0][1][3] is the X-coordinate of subap number 3 of wfs number 0 at altitude 3

  //Computes  u and v
  sub_pos_gpu_gb(tomo_gpu, tomo);
  
  cudaStreamSynchronize(tomo_gpu->matcov_stream);
  if (sspSizeL) free(sspSizeL);
}
//---------------------------------------------------------------------------------
void update_tomo_sys_gpu_gb(struct gtomo_struct *tomo_gpu,struct tomo_struct tomo) 
{

  cudaError_t e;

  long ioff[tomo.Nw];
  ioff[0] = 0;
  for (int i=1;i<tomo.Nw;i++) ioff[i] = ioff[i-1] + tomo.Nsubap[i-1];
  e = cudaMemcpyAsync(tomo_gpu->ioff_d, ioff, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu ioff_d");

  e = cudaMemcpyAsync(tomo_gpu->alphaX_d, tomo.alphaX, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu alphaX_d");
  e = cudaMemcpyAsync(tomo_gpu->alphaY_d, tomo.alphaY, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu alphaY_d");

  e = cudaMemcpyAsync(tomo_gpu->GsAlt_d, tomo.GsAlt, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu GsAlt_d");

  e = cudaMemcpyAsync(tomo_gpu->Nssp_d, tomo.Nssp, tomo.Nw*sizeof(long), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu Nssp_d");

  e = cudaMemcpyAsync(tomo_gpu->diamPup_d, tomo.diamPup, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu diamPup_d");

  e = cudaMemcpyAsync(tomo_gpu->XPup_d, tomo.XPup, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu XPup_d");
  e = cudaMemcpyAsync(tomo_gpu->YPup_d, tomo.YPup, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu YPup_d");
  e = cudaMemcpyAsync(tomo_gpu->thetaML_d, tomo.thetaML, tomo.Nw*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu thetaML_d");

  e = cudaMemcpyAsync(tomo_gpu->X_d, tomo.X, tomo.Nx*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu X_d");
  e = cudaMemcpyAsync(tomo_gpu->Y_d, tomo.Y, tomo.Nx*sizeof(double), cudaMemcpyHostToDevice, tomo_gpu->matcov_stream);
  process_err(e, "copy gpu Y_d");
  
  //cudaStreamSynchronize(tomo_gpu->matcov_stream);
}


//============================================================================================
//============================= MATCOV 3 =====================================================
//============================================================================================
//extern "C"
void matcov_gpu_3(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu)
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
	
  //cudaError_t e;
  
  char uplo, copy; 
  
  uplo = 'f';	// full generation is enabled by default
  copy = 'c';
  
  int type_mat = tomo.part;
  
  if(type_mat == 1) // Caa matrix
  {
  	// check if a square diagonal tile is generated then we set uplo to 'l' or 'u'
  	// and then enable the copy
  	// This also applies if the entire matrix will be generated
  	// otherwise (off diagonal tile or non square submatrix) - full generation is assumed
  	if((xoffset == yoffset) && (nrows == ncols))	// if sqaure & diagonal
  	{
  		uplo = 'l'; 
  		if(type_mat == 1)copy = 'c';
  	}
  	else	// full generation, copy is ignored
  	{
  		uplo = 'f';
  	}
  }
  else if(type_mat == 2 || type_mat == 3) // Cmaa matrix
  {
  	uplo = 'f';		// full generation, copy is ignored
  }
  else
  {
  	printf("ERROR: unrecognized type_mat %d \n", type_mat); exit(1);
  }
  
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
  //const long Nw = tomo.Nw;
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);

  //int size = tomo.Nslopes - 2 * tomo.Nsubap[tomo.Nw-1];

  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = tomo.Nsubap[0];
  
  // generate a full matrix
  matcov_kernel_3<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(uplo, copy, data, nrows, ncols, xoffset, yoffset, lda,
					   convert,tomo_gpu->sspSizeL_d,tomo_gpu->Nssp_d,tomo_gpu->u_d,tomo_gpu->v_d,
					   pasDPHI,tomo_gpu->tabDPHI_d,tomo_gpu->indexL0_d,tomo_gpu->cn2_d,
					   Ndphi,tomo.Nw,tomo.Nlayer,Nsubap,type_mat,tomo.DiamTel);
  
  //if (type_mat == 1)
  //  matcov_kernel_copy<<<dimGrid, dimBlock>>>(data, nrows, ncols, xoffset, yoffset, lda);
  
  //cudaStreamSynchronize(tomo_gpu->matcov_stream);
}



//============================================================================================
//=========================== MATCOV 4 (NOISE) KERNELS/FUNCTION ==============================
//============================================================================================
__device__ double compute_element_4(int ipos, int jpos, double convert, double *sspSizeL, long *Nssp, double *u, double *v, 
					double pasDPHI, double *tabDPHI, long *indexL0, double *cn2, int Ndphi, int Nw, int Nlayer, 
					int Nsubap, double *alphaX, double *alphaY, double lgs_cst, double noise_var, double spotWidth,
					double dH_lgs, double alt_lgs, int type_mat, int nlgs, double teldiam)
{
	/* *** Covariance matrix per-element generation ***
	*   Arguments
	*   =========
	*	ipos:		Integer: global x-coordinate of the element w.r.t. the entire matrix
	*	jpos:		Integer: global y-coordinate of the element w.r.t. the entire matrix
	*/
	
	// for now return a dummy value
  
  const double lambda2 = 0.00026942094446267851;
  //WFS m
  int m = ipos / (2 * Nsubap);
  if (type_mat == 3) m = Nw-1;
  //WFS n
  int n = jpos / (2 * Nsubap);
  if (type_mat == 2) n = Nw-1;
  //subap i
  int i = ipos % (2 * Nsubap); 
  //subap j
  int j = jpos % (2 * Nsubap);
  //xy i
  int xy_i;
  //xy j
  int xy_j;
  if (i>=Nsubap) {
    i-= Nsubap;
    xy_i = 1;
  } else xy_i = 0;
  if (j>=Nsubap) {
    j-= Nsubap;
    xy_j = 1;
  } else xy_j = 0;

  const double sspSizem = teldiam / Nssp[m];
  const double sspSizen = teldiam / Nssp[n];
  
  const double kk = lambda2 / (sspSizem * sspSizen);
    
  int type = xy_i * 2 + xy_j;

  //Layer l
  double covar = 0.0;
  #pragma unroll
  for (int l = 0; l < Nlayer; l++) 
  {
    double sspSizeml = sspSizeL[m * Nlayer + l];
    double sspSizenl = sspSizeL[n * Nlayer + l];
    //test if the altitude layers is not higher than the LGS altitude
    if ((sspSizeml > 0) && (sspSizenl > 0)) 
    {
      int pos1 = m + i * Nw + l * Nw * Nsubap;
      int pos2 = n + j * Nw + l * Nw * Nsubap;
      //if(threadIdx.x == 6 && threadIdx.y == 0 && blockIdx.x == 6 && blockIdx.y == 1)
      //if((pos1 >= 6840) || (pos2 >= 6839))
      //{
      //	printf("================ pos1 = %d, pos2 = %d \n", pos1, pos2);
      //}
      //(6,0,0) in block (0,2,0);
      double du =  u[pos1] - u[pos2];	      
      double dv =  v[pos1] - v[pos2];
      
      double s1 = sspSizeml * 0.5;
      double s2 = sspSizenl * 0.5;
      
      double ac = s1 - s2;
      double ad = s1 + s2;
      double bc = -ad;   // initially -s1-s2;
      double bd = -ac;   // initially -s1+s2;

      if (type == 0) covar += 0.5 * pasDPHI * cov_XX_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else if (type == 3) covar += 0.5 * pasDPHI * cov_YY_gpu_gb(du,dv,ac,ad,bc,bd,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l];
      else //if ((type == 1) || (type == 2)) 
      {
      	double s0 = sqrt(s1 * s1 + s2 * s2); //half size of the subaperture equivalent to a convolution by s1 and s2
      	double dd = (s1 > s2) ? 1. - s2 / s1 : 1. - s1 / s2; // Nono's style ....
      	covar += 0.25 * pasDPHI * cov_XY_gpu_gb(du,dv,s0,tabDPHI,indexL0[l],convert,Ndphi) * kk * cn2[l] * (1. - dd * dd);
      }
    }
  }
  // adding noise
  if (m == n) {
    if (m < nlgs) {
      if (i == j) {
	// lgs case
	const int pos1 = m + i * Nw;
	double x = u[pos1];	      
	double y = v[pos1];
	double xwfs = alphaX[m] * 206265;	      
	double ywfs = alphaY[m] * 206265;
	double lltx = 0;	      
	double llty = 0;
	const double lltnorm = sqrtf(xwfs*xwfs + ywfs*ywfs);
	if (lltnorm != 0) {
	  lltx = xwfs / lltnorm * teldiam / 2.0;
	  llty = ywfs / lltnorm * teldiam / 2.0;
	}
	x -= lltx;
	y -= llty;
        x  = 206265. * dH_lgs * x / alt_lgs / alt_lgs;   // extension at Fwhm, in arcsec
        y  = 206265. * dH_lgs * y / alt_lgs / alt_lgs;   // extension at Fwhm, in arcsec
        double lgsExt = sqrtf(x * x + y * y);   // lengh of the extension
        double lgsTheta = x != 0 ? atanf( y / x) : 0.0;   // angle of extension
        double totalExt = sqrtf( lgsExt *  lgsExt + spotWidth * spotWidth); 
	// lengh of the extension including seeing, laser size, ...
	double ratio = totalExt / spotWidth;
        double noiseLongAxis = noise_var * ratio * ratio;
	if (type == 0) covar += noiseLongAxis * cosf(lgsTheta) * cosf(lgsTheta) + 
			 noise_var * sinf(lgsTheta) * sinf(lgsTheta);
	else if (type == 3) covar += noiseLongAxis * sinf(lgsTheta) * sinf(lgsTheta) + 
			      noise_var * cosf(lgsTheta) * cosf(lgsTheta);
	else covar += (noiseLongAxis-noise_var) * sinf(lgsTheta) * cosf(lgsTheta);
      }
      if ((type == 0) || (type == 3))
	covar += lgs_cst;
    } else {
    // ngs case
      if (i==j) {
	if ((type == 0) || (type == 3)) {
	  covar += noise_var;
	}
      }
    }
  }

  return (double)covar; 
}

//------------------------------------------------------------------------------------------
__global__ void matcov_kernel_4(char uplo, char copy, double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
				       double convert, double *sspSizeL, long *Nssp, double *u, double *v, 
					double pasDPHI, double *tabDPHI, long *indexL0, double *cn2, int Ndphi, int Nw, int Nlayer, 
					int Nsubap, double *alphaX, double *alphaY, double lgs_cst, double noise_var, double spotWidth,
					double dH_lgs, double alt_lgs, int type_mat, int nlgs, double teldiam)
{
  /* *** covariance matrix generation kernel ***
   *	The kernel generates the element values in a given matrix/submatrix
   *   The generation function can be any function, as long as each element 
   *   can be computed both individually and independently
   *
   *	see argument description in the kernel driver
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
	
  double value;
  if(uplo == 'l')
  {
  	if(gy <= gx)
  	{	
  		value = compute_element_4(	gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, 
							Nsubap, alphaX, alphaY, lgs_cst, noise_var, spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
  		data[ly * lda + lx] = value;
  		if(copy == 'c') data[lx * lda + ly] = value;
  	}
  }
  else if (uplo == 'u') // upper
  {
  	if(gx <= gy)
  	{
  		value = compute_element_4(	gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, 
							Nsubap, alphaX, alphaY, lgs_cst, noise_var, spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
  		data[ly * lda + lx] = value;
  		if(copy == 'c') data[lx * lda + ly] = value;
  	}
  }
  else	// uplo = 'f' full generation
  {
  	value = compute_element_4(	gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, 
							Nsubap, alphaX, alphaY, lgs_cst, noise_var, spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
  	data[ly * lda + lx] = value;
  }
  
  //if ((type_mat == 3) || (gx <= gy)) {
  //  // call the generation function
  //  data[0] = compute_element_4(gx, gy, convert, sspSizeL, Nssp, u, v, pasDPHI, tabDPHI, indexL0, cn2, Ndphi, Nw, Nlayer, 
	//				Nsubap, alphaX, alphaY, lgs_cst, noise_var, spotWidth, dH_lgs, alt_lgs, type_mat, nlgs, teldiam);
    //printf("gx = %d, gy = %d ----- %.2f \n", gx, gy, data[0]);
  //} 
}

//============================================================================================
//============================= MATCOV 4 (NOISE) =============================================
//============================================================================================
void matcov_gpu_4(double* data, int nrows, int ncols, int xoffset, int yoffset, int lda, struct tomo_struct tomo, struct gtomo_struct *tomo_gpu)
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
	
  //cudaError_t e;
  char uplo, copy; 
  
  uplo = 'f';	// full generation is enabled by default
  copy = 'c';
  
  int type_mat = tomo.part;
  
  if(type_mat == 1) // Caa matrix
  {
  	// check if a square diagonal tile is generated then we set uplo to 'l' or 'u'
  	// and then enable the copy
  	// This also applies if the entire matrix will be generated
  	// otherwise (off diagonal tile or non square submatrix) - full generation is assumed
  	if((xoffset == yoffset) && (nrows == ncols))	// if sqaure & diagonal
  	{
  		uplo = 'l'; 
  		copy = 'c';
  	}
  	else	// full generation, copy is ignored
  	{
  		uplo = 'f';
  	}
  }
  //else if(type_mat == 2) //
  else if(type_mat == 2 || type_mat == 3) // Cmaa matrix
  {
  	uplo = 'f';		// full generation, copy is ignored
  }
  else
  {
  	printf("ERROR: unrecognized type_mat %d \n", type_mat); exit(1);
  }
  // %%%%%%% Pre-computation of DPHI %%%%%%%%%%
  //Computes an array of DPHI (tabDPHI) for an array of subaperture distance rr for each DIFFERENT L0
  const double crmax = tomo.rmax;
  const double pasDPHI = 1./tomo.pasDPHI; //inverse du pas de rr
  const long Ndphi = floor(crmax*pasDPHI)+1;
  const double convert = (double)(Ndphi-1)/(crmax+1./pasDPHI);

  int nbx = nrows / matcov_thread_x + (nrows%matcov_thread_x != 0);
  int nby = ncols / matcov_thread_y + (ncols%matcov_thread_y != 0);

  dim3 dimBlock(matcov_thread_x, matcov_thread_y);
  dim3 dimGrid(nbx, nby);
  const long Nsubap = tomo.Nsubap[0];
  
  matcov_kernel_4<<<dimGrid, dimBlock, 0, tomo_gpu->matcov_stream>>>(uplo, copy, data, nrows, ncols, xoffset, yoffset, lda, convert, tomo_gpu->sspSizeL_d, 
						tomo_gpu->Nssp_d, tomo_gpu->u_d, tomo_gpu->v_d, pasDPHI, tomo_gpu->tabDPHI_d, 
						tomo_gpu->indexL0_d, tomo_gpu->cn2_d, Ndphi, tomo.Nw, tomo.Nlayer, 
						Nsubap, tomo_gpu->alphaX_d, tomo_gpu->alphaY_d, tomo.lgs_cst, tomo.noise_var, 
						tomo.spot_width, tomo.lgs_depth, tomo.lgs_alt, type_mat, tomo.nlgs, tomo.DiamTel);

	//cudaStreamSynchronize(tomo_gpu->matcov_stream);
	
  //if (type_mat == 1)
  // matcov_kernel_copy<<<dimGrid, dimBlock>>>(data, nrows, ncols, xoffset, yoffset, lda);
}
