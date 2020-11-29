/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "matcov_kernels.hpp"
#include <math.h>
#include "moao_defs.h"


/*! \brief Computation of the function

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

 */
 template<typename T>
T macdo_x56(T x, int k)
{
  const double a = 5. / 6.;
  const T x2a = pow(x, (T)2. * a);
  const T x22 = x * x / 4.;
  T x2n;               // x^2.a, etc
  T s = 0.0;
  int n;

  x2n = 0.5;                           // init (1/2) * x^0

  s = MOAO_CST::Gma[0] * x2a;
  s *= x2n;

  // prepare recurrence iteration for next step
  x2n *= x22;    // x^n

  for (n = 1; n <= 10; n++) {
    s += (MOAO_CST::Gma[n] * x2a + MOAO_CST::Ga[n]) * x2n;
    // prepare recurrence iteration for next step
    x2n *= x22;    // x^n
  }

  return s;
}

template<typename T>
T asymp_macdo(T x)
/* DOCUMENT asymp_macdo(x)

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
  T res;
  T x_1;

  x_1 = 1. / x;
  res = MOAO_CST::k2
      - MOAO_CST::k3 * exp(-x) * pow(x, (T)(1 / 3.))
          * (1.0 + x_1 * (MOAO_CST::a1 + x_1 * (MOAO_CST::a2 + x_1 * MOAO_CST::a3)));
  return res;
}

template<typename T>
T rodconan(T r, T L0, int k)
/* DOCUMENT rodconan(r,L0,k=)
 The phase structure function is computed from the expression
 Dphi(r) = k1  * L0^(5./3) * (k2 - (2.pi.r/L0)^5/6 K_{5/6}(2.pi.r/L0))

 For small r, the expression is computed from a development of
 K_5/6 near 0. The value of k2 is not used, as this same value
 appears in the series and cancels with k2.
 For large r, the expression is taken from an asymptotic form.

 SEE ALSO:
 */
{
  T res = 0;

  const T dprf0 = (2 * MOAO_CST::pi / L0) * r;

  if (dprf0 > MOAO_CST::Xlim) { 
    res = asymp_macdo(dprf0);
  } else {
    res = -macdo_x56(dprf0, k);
  }
  res *= MOAO_CST::k1 * pow(L0, (T)5. / (T)3.0);
  return res;
}

//============================================================================================
/*!
 * dphi = DPHI(x,y,indexL0,rr,tabDPHI,convert) * r0^(-5./3)
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
template<typename T>
T DPHI_gb(T x, T y, T indexL0)//, T *tabDPHI, T convert, int Ndphi)
{
  T r = sqrt(x * x + y * y);
  
  return rodconan(r, indexL0, 10);
}

//------------------------------------------------------------------------------------
/*!
 *  Compute the XX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
template<typename T>
T cov_XX_gb(T du, T dv, T ac, T ad, T bc, T bd, T indexL0)

{
  return -DPHI_gb(du + ac, dv, indexL0)
    + DPHI_gb(du + ad, dv, indexL0)
    + DPHI_gb(du + bc, dv, indexL0)
    - DPHI_gb(du + bd, dv, indexL0);;
}

//------------------------------------------------------------------------------------
/*!
 *  Compute the YY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
template<typename T>
T cov_YY_gb(T du, T dv, T ac, T ad, T bc, T bd, T indexL0)
{ 
  return  -DPHI_gb(du, dv + ac, indexL0)
    + DPHI_gb(du, dv + ad, indexL0)
    + DPHI_gb(du, dv + bc, indexL0)
    - DPHI_gb(du, dv + bd, indexL0);;
}

//------------------------------------------------------------------------------------
/*!
 * Compute the XY-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
template<typename T>
T cov_XY_gb(T du, T dv, T s1, T s2, T indexL0)
{
  return -DPHI_gb(du + s1 , dv - s2, indexL0)
        + DPHI_gb(du + s1, dv + s2, indexL0)
        + DPHI_gb(du - s1, dv - s2, indexL0)
        - DPHI_gb(du - s1, dv + s2, indexL0);
}

//------------------------------------------------------------------------------------
/*!
 * Compute the YX-covariance with the distance sqrt(du2+dv2). DPHI is precomputed on tabDPHI.
 */
template<typename T>
T cov_YX_gb(T du, T dv, T s1, T s2, T indexL0)
{
  return -DPHI_gb(du - s2, dv + s1, indexL0)
        + DPHI_gb(du + s2, dv + s1, indexL0)
        + DPHI_gb(du - s2, dv - s1, indexL0)
        - DPHI_gb(du + s2, dv - s1, indexL0);
}




//------------------------------------------------------------------------------------
/*! Covariance matrix per-element generation ***
*   Arguments
*   =========
*       ipos:           Integer: global x-coordinate of the element w.r.t. the entire matrix
*       jpos:           Integer: global y-coordinate of the element w.r.t. the entire matrix
*/
template<typename T>
T compute_element_tiled_4(
    int ipos, int jpos, T *sspSizeL, long *Nssp, T *u, T *v,
    T *indexL0, T *cn2, int Nw, int Nlayer,
    long * Nsubap_wfs, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx,
    T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T teldiam)
{
  
  const T lambda2 = 0.00026942094446267851;
  //WFS m
  //int m = ipos / (2 * Nsubap);
  long Nsubapx = Nsubap_wfs[0];
  int m = 0;
  if (type_mat == 3){
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
  //int n = jpos / (2 * Nsubap);
  long Nsubapy = Nsubap_wfs[0];
  int n = 0;
  if (type_mat == 2){
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
  int l;
  for (l = 0; l < Nlayer; l++) 
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
          covar += 0.5 * cov_XX_gb(du,dv,ac,ad,bc,bd,indexL0[l]) * kk * cn2[l];
      }
      else if (type == 3){
          covar += 0.5 * cov_YY_gb(du,dv,ac,ad,bc,bd,indexL0[l]) * kk * cn2[l];
      }
      else if(type==1){
          covar += 0.5 * cov_XY_gb(du,dv,s1,s2,indexL0[l]) * kk * cn2[l];
      }
      else if(type==2){
          covar += 0.5 * cov_YX_gb(du,dv,s1,s2,indexL0[l]) * kk * cn2[l];
      }
    }
  }

  // adding noise
  if (m == n) {
    if (m < nlgs) {
      if (i == j) {
        // lgs case
        const int pos1 = i + Nsubapx;
        if (type == 0) covar += noiseLGSxx[pos1];
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

  return (T)covar; 
}
template float compute_element_tiled_4(
    int ipos, int jpos, float *sspSizeL, long *Nssp, float *u, float *v,
    float *indexL0, float *cn2, int Nw, int Nlayer,
    long * Nsubap_wfs, long Nx, float lgs_cst, float *noiseNGS, float *noiseLGSxx,
    float *noiseLGSyy, float *noiseLGSxy, int type_mat, int nlgs, float teldiam);
template double compute_element_tiled_4(
    int ipos, int jpos, double *sspSizeL, long *Nssp, double *u, double *v,
    double *indexL0, double *cn2, int Nw, int Nlayer,
    long * Nsubap_wfs, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx,
    double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double teldiam);

//============================================================================================
/*!
 * covariance matrix generation kernel ***
 *    The kernel generates the element values in a given matrix/submatrix
 *   The generation function can be any function, as long as each element
 *   can be computed both individually and independently
 *
 *    see argument description in the kernel driver
 */
 template<typename T>
void matcov_kernel_4(
  char uplo, char copy, T* data, int xoffset, int yoffset, int lda,
  T *sspSizeL, long *Nssp, T *u, T *v,
  T *indexL0, T *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, T lgs_cst, T *noiseNGS, T *noiseLGSxx, 
  T *noiseLGSyy, T *noiseLGSxy, int type_mat, int nlgs, T teldiam, int lx, int ly)
{
        
  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;
  T value;

  if(uplo == 'l')
  {
    if(gy <= gx)
    {
      value = compute_element_tiled_4(
                  gx, gy, sspSizeL, Nssp, u, v, indexL0, cn2,
                  Nw, Nlayer, Nsubap, Nx, lgs_cst,noiseNGS, noiseLGSxx,
                  noiseLGSyy, noiseLGSxy, type_mat, nlgs, teldiam);
      data[ly * lda + lx] = value;
      if(copy == 'c') data[lx * lda + ly] = value;
    }
  }
  else if (uplo == 'u') // upper
  {
    if(gx <= gy)
    {
      value = compute_element_tiled_4(
                  gx, gy, sspSizeL, Nssp, u, v, indexL0, cn2,
                  Nw, Nlayer, Nsubap, Nx, lgs_cst, noiseNGS, noiseLGSxx,
                  noiseLGSyy, noiseLGSxy, type_mat, nlgs, teldiam);
      data[ly * lda + lx] = value;
      if(copy == 'c') data[lx * lda + ly] = value;
    }
  }
  else  // uplo = 'f' full generation
  {
    value = compute_element_tiled_4(
                  gx, gy, sspSizeL, Nssp, u, v, indexL0, cn2,
                  Nw, Nlayer, Nsubap, Nx, lgs_cst, noiseNGS, noiseLGSxx,
                  noiseLGSyy, noiseLGSxy, type_mat, nlgs, teldiam);
    data[ly * lda + lx] = value;
  }
}
template void matcov_kernel_4(
  char uplo, char copy, float* data, int xoffset, int yoffset, int lda,
  float *sspSizeL, long *Nssp, float *u, float *v,
  float *indexL0, float *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, float lgs_cst, float *noiseNGS, float *noiseLGSxx, 
  float *noiseLGSyy, float *noiseLGSxy, int type_mat, int nlgs, float teldiam, int lx, int ly);
template void matcov_kernel_4(
  char uplo, char copy, double* data, int xoffset, int yoffset, int lda,
  double *sspSizeL, long *Nssp, double *u, double *v,
  double *indexL0, double *cn2, int Nw, int Nlayer,
  long *Nsubap, long Nx, double lgs_cst, double *noiseNGS, double *noiseLGSxx, 
  double *noiseLGSyy, double *noiseLGSxy, int type_mat, int nlgs, double teldiam, int lx, int ly);

//============================================================================================
/*!
 * ** Covariance matrix per-element generation ***
 *   Arguments
 *   =========
 *       ipos:           Integer: global x-coordinate of the element w.r.t. the entire matrix
 *       jpos:           Integer: global y-coordinate of the element w.r.t. the entire matrix
 */
 template<typename T>
T compute_element_ts_tile(
  int ipos, int jpos,T *X, T *Y,
  long *Nssp, T *L0diff, T *cn2, 
  int Nw, int Nlayer, int Nsubap, T teldiam)
{
  //(lambda/(2*pi)*rasc)**2  with lambda=0.5e-6 m
  //TODO generalize lambda
  T lambda2 = 0.00026942094446267851;
  //WFS Nw-1
   //subap i
  int i = ipos < Nsubap ? ipos : ipos - Nsubap;
  //subap j
  int j = jpos < Nsubap ? jpos : jpos - Nsubap;
  //xy i
  int xy_i = ipos < Nsubap ? 0 : 1;
  //xy j
  int xy_j = jpos < Nsubap ? 0 : 1;
  
  T sspSize = teldiam / Nssp[Nw-1];
  
  T kk = lambda2 / (sspSize * sspSize);
    
  int type = xy_i * 2 + xy_j;

  T s = sspSize * 0.5;
  
  T ac = 0.0;
  T ad = 2.0 * s;
  T bc = -ad;   
  T bd = 0.0;   

    //TODO: valable uniquement si Nsubap constant
  T du = X[(Nsubap*(Nw-1)+i)] - X[(Nsubap*(Nw-1)+j)];            
  T dv = Y[(Nsubap*(Nw-1)+i)] - Y[(Nsubap*(Nw-1)+j)];
  
  //Layer l
  T covar = 0.0;
  //#pragma unroll
  int l;
  for (l = 0; l < Nlayer; l++) 
  {
     //test if the altitude layers is not higher than the LGS altitude
    if (sspSize > 0) 
    {
      if (type == 0) {
          covar += 0.5 * cov_XX_gb(du,dv,ac,ad,bc,bd,L0diff[l]) * kk * cn2[l];
      }
      else if (type == 3){
          covar += 0.5 * cov_YY_gb(du,dv,ac,ad,bc,bd,L0diff[l]) * kk * cn2[l];
      }
      else if(type==1){
          covar += 0.5 * cov_XY_gb(du,dv,s,s,L0diff[l]) * kk * cn2[l];
      }
      else if(type==2){
          covar += 0.5 * cov_YX_gb(du,dv,s,s,L0diff[l]) * kk * cn2[l];
      }

    }
  }
  return (T)covar; 
}


//--------------------------------------------------------------------------------------------
/*!
 *** covariance matrix generation kernel ***
 *    The kernel generates the element values in a given matrix/submatrix
 *   The generation function can be any function, as long as each element
 *   can be computed both individually and independently
 *
 *    see argument description in the kernel driver
 */
 template<typename T>
void matcov_ts_kernel_tile(
  T* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  T *X, T *Y, long *Nssp,
  T *indexL0, T *cn2, int Nw, int Nlayer, int Nsubap, T teldiam,
  int lx, int ly)
{
        
  // global coordinates of the elemnt w.r.t. the entire matrix
  int gx = lx + xoffset;
  int gy = ly + yoffset;

  // out-of-bound threads should terminate
  if( (lx >= nrows) || (ly >= ncols) ) return;
        
  // Advance the data pointer accordingly
  data += ly * lda + lx;
        
  // call the generation function
  data[0] = compute_element_ts_tile(gx, gy,X, Y,Nssp, indexL0,cn2,Nw,Nlayer,Nsubap,teldiam);
}
template void matcov_ts_kernel_tile(
  float* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  float *X, float *Y, long *Nssp, float *indexL0, float *cn2, int Nw, 
  int Nlayer, int Nsubap, float teldiam, int lx, int ly);
template void matcov_ts_kernel_tile(
  double* data, int nrows, int ncols, int xoffset, int yoffset, int lda,
  double *X, double *Y, long *Nssp, double *indexL0, double *cn2, int Nw, 
  int Nlayer, int Nsubap, double teldiam, int lx, int ly);

