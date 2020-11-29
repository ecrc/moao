#include "matcov.hpp"

#include "codelet_matcov.hpp"
#include "chameleon_templates.hpp"
#include <starpu.h>
#include <morse.h>
#include "context.h"


/*! matcov tile driver
*  Arguments
*  ==========
*  data         double pointer: A pointer to the matrix/submatrix to be generated. It
*                       should always point to the first element in a matrix/submatrix
*
*  nrows        integer: The number of rows of the matrix/submatrix to be generated
*
*  ncols        integer: The number of columns of the matrix/submatrix to be generated
*
*  xoffset      integer: The x-offset of the submatrix, must be zero if the entire matrix
*                       is generated. Its the x-coordinate of the first element in the matrix/submatrix
*
*  yoffset  integer: The y-offset of the submatrix, must be zero if the entire matrix
*                       is generated. Its the y-coordinate of the first element in the matrix/submatrix
*
*  lda          integer: The leading dimension of the matrix/submatrix
*
*  rowMajor     integer: specify if the tile is generated in row major or column major layout
*                        0:column major, else row major (default=1)
*/
template<typename T>
 int MOAO_CHAMELEON::matcov_tile(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, Tomo_struct<T> *tomo, int type_mat)
{
MORSE_desc_t *descSspSizeL=NULL;
MORSE_desc_t *descNssp=NULL;
MORSE_desc_t *descU=NULL;
MORSE_desc_t *descV=NULL;
MORSE_desc_t *descX=NULL;
MORSE_desc_t *descY=NULL;
MORSE_desc_t *descL0diff=NULL;
MORSE_desc_t *descCn2=NULL;
MORSE_desc_t *descNsubap=NULL;
MORSE_desc_t *descAlphaX=NULL;
MORSE_desc_t *descAlphaY=NULL;
MORSE_desc_t *descNoiseNGS;
MORSE_desc_t *descNoiseLGSxx;
MORSE_desc_t *descNoiseLGSyy;
MORSE_desc_t *descNoiseLGSxy;


int Nw=(int)tomo->sys.nW, Nl=(int)tomo->atm.nLayer, Nx=(int)tomo->Nx;
  int i;
  int Nsubap=0;
  for(i=0;i<tomo->sys.nLgs;i++){
  Nsubap+=tomo->Nsubap[i];
  }
  if(!Nsubap)
    Nsubap=1;

MORSE_Desc_Create(&(descSspSizeL)   , tomo->sspSizeL.data(),    MORSE<T>::real , Nw*Nl  , 1, Nw*Nl , Nw*Nl , 1, 0, 0, Nw*Nl , 1,1,1);
MORSE_Desc_Create(&(descNssp)       , tomo->Nssp.data(),        MORSE<T>::int64, Nw     , 1, Nw    , Nw    , 1, 0, 0, Nw    , 1,1,1);
MORSE_Desc_Create(&(descU)          , tomo->u.data() ,          MORSE<T>::real , Nl*Nx  , 1, Nl*Nx , Nl*Nx , 1, 0, 0, Nl*Nx , 1,1,1);
MORSE_Desc_Create(&(descV)          , tomo->v.data() ,          MORSE<T>::real , Nl*Nx  , 1, Nl*Nx , Nl*Nx , 1, 0, 0, Nl*Nx , 1,1,1);
MORSE_Desc_Create(&(descX)          , tomo->X.data(),           MORSE<T>::real , Nx     , 1, Nx    , Nx    , 1, 0, 0, Nx    , 1,1,1);
MORSE_Desc_Create(&(descY)          , tomo->Y.data(),           MORSE<T>::real , Nx     , 1, Nx    , Nx    , 1, 0, 0, Nx    , 1,1,1);
MORSE_Desc_Create(&(descL0diff)     , tomo->L0diff.data(),      MORSE<T>::real , Nl     , 1, Nl    , Nl    , 1, 0, 0, Nl    , 1,1,1);
MORSE_Desc_Create(&(descCn2)        , tomo->atm.cn2.data(),     MORSE<T>::real , Nl     , 1, Nl    , Nl    , 1, 0, 0, Nl    , 1,1,1);
MORSE_Desc_Create(&(descNsubap)     , tomo->Nsubap.data(),      MORSE<T>::int64, Nw     , 1, Nw    , Nw    , 1, 0, 0, Nw    , 1,1,1);
MORSE_Desc_Create(&(descAlphaX)     , tomo->alphaX.data(),      MORSE<T>::real , Nw     , 1, Nw    , Nw    , 1, 0, 0, Nw    , 1,1,1);
MORSE_Desc_Create(&(descAlphaY)     , tomo->alphaY.data(),      MORSE<T>::real , Nw     , 1, Nw    , Nw    , 1, 0, 0, Nw    , 1,1,1);
MORSE_Desc_Create(&(descNoiseNGS)   , tomo->noiseNGS.data(),    MORSE<T>::real , Nw     , 1, Nw    , Nw    , 1, 0, 0, Nw    , 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSxx) , tomo->noiseLGSxx.data(),  MORSE<T>::real , Nsubap , 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSyy) , tomo->noiseLGSyy.data(),  MORSE<T>::real , Nsubap , 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);
MORSE_Desc_Create(&(descNoiseLGSxy) , tomo->noiseLGSxy.data(),  MORSE<T>::real , Nsubap , 1, Nsubap, Nsubap, 1, 0, 0, Nsubap, 1,1,1);

    MORSE_context_t *morse;
    MORSE_option_t options;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return sequence->status;
    RUNTIME_options_init(&options, morse, sequence, request);

    int m,n,ldam;
    int tempmm,tempnn;
    int ret=starpu_init(NULL);

    for(m=0;m<descA->mt;m++){
        tempmm = m == descA->mt-1 ? descA->m-m*descA->mb : descA->mb;
        ldam=descA->get_blkldd(descA, m);
        for(n=0;n<descA->nt;n++){
            tempnn = n == descA->nt-1 ? descA->n-n*descA->nb : descA->nb;
            MOAO_STARPU::TASK_matcov(descA,tempmm,tempnn,m,n,ldam,*tomo,
            descSspSizeL,
            descNssp,
            descU,
            descV,
            descX,
            descY,
            descL0diff,
            descCn2,
            descNsubap,
            descNoiseNGS,
            descNoiseLGSxx,
            descNoiseLGSyy,
            descNoiseLGSxy,
            type_mat);
        }
    }

    RUNTIME_options_finalize(&options, morse);
    MORSE_Desc_Flush(descA          ,sequence);


MORSE_Desc_Destroy(&descSspSizeL);
MORSE_Desc_Destroy(&descNssp);
MORSE_Desc_Destroy(&descU);
MORSE_Desc_Destroy(&descV);
MORSE_Desc_Destroy(&descX);
MORSE_Desc_Destroy(&descY);
MORSE_Desc_Destroy(&descL0diff);
MORSE_Desc_Destroy(&descCn2);
MORSE_Desc_Destroy(&descNsubap);
MORSE_Desc_Destroy(&descAlphaX);
MORSE_Desc_Destroy(&descAlphaY);
MORSE_Desc_Destroy(&descNoiseNGS);
MORSE_Desc_Destroy(&descNoiseLGSxx);
MORSE_Desc_Destroy(&descNoiseLGSyy);
MORSE_Desc_Destroy(&descNoiseLGSxy);

    return 0;
}

template int MOAO_CHAMELEON::matcov_tile(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, Tomo_struct<float> *tomo, int type_mat);
template int MOAO_CHAMELEON::matcov_tile(MORSE_desc_t *descA, MORSE_sequence_t *sequence, MORSE_request_t *request, Tomo_struct<double> *tomo, int type_mat);
