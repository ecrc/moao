
cimport numpy as np
import numpy as np

include "par.pxi"

include "tomo.pyx"
include "intersample.pyx"
import commande as cmd


IF USE_SINGLE:
    cython_real=np.float32
ELSE:
    cython_real=np.float64

#
#pipeline routines
#
import matplotlib.pyplot as pl
pl.ion()

#def py_reconstructor(np.ndarray[ndim=2,dtype=cython_real_t] cmm, np.ndarray[ndim=2,dtype=cython_real_t] ctm,int M, int N,int ngal=1):
#
#    cdef np.ndarray[ndim=2,dtype=cython_real_t] tmpCmm,tmpCtm
#
#    tmpCmm=np.asfortranarray(cmm) # since we use COL_MAJOR and lower part, transpose replace upper by lower
#    tmpCtm=np.asfortranarray(ctm)
#
#    cdef real_t *cmm_data=<real_t*>tmpCmm.data
#    cdef real_t *ctm_data=<real_t*>tmpCtm.data
#
#    reconstructor(M,N,cmm_data,cmm.shape[1],ctm_data,ctm.shape[0],ngal)
#
#    return tmpCtm
#
#
#def py_computeCee_Cvv(np.ndarray[ndim=2,dtype=cython_real_t] Cmm, np.ndarray[ndim=2,dtype=cython_real_t] Cpp,
#                      np.ndarray[ndim=2,dtype=cython_real_t] Cpm, np.ndarray[ndim=2,dtype=cython_real_t] R,
#                      np.ndarray[ndim=2,dtype=cython_real_t] Dx):
#
##   int nMeasTS: number of measurements of the truth sensor
##   int nMeas   :number of measurements
##   real_t *Cmm :covariance matrix between the sensors                              Cmm(nMeas,nMeas)
##   real_t *Cpp :covariance matrix of the truth sensor                              Cpp(nMeasTS,nMeasTS)
##   real_t *Cpm :covariance matrix between the truth sensors and the other sensors  Cpm(nMeas,nMeasTS)
##   real_t *R   :tomographic reconstructor                                          R(nMeas,nMeasTS))
##   real_t *Dx  :                                                                   Dx(nact,nMeasTS)
#
#    #check dimension of all input: TODO
#    cdef int nMeas   =Cmm.shape[0]
#    cdef int nMeasTS =Cpp.shape[0]
#    cdef int nact    =Dx.shape[0]
#
#    cdef np.ndarray[ndim=2,dtype=cython_real_t] cmm,cpp,cpm,r,dx
#    cpp=np.copy(np.asfortranarray(Cpp))
#    cpm=np.copy(np.asfortranarray(Cpm))
#    r  =np.copy(np.asfortranarray(R  ))
#    dx =np.copy(np.asfortranarray(Dx ))
#
#    cdef np.ndarray[ndim=2,dtype=cython_real_t] Cee=np.asfortranarray(np.zeros((nMeasTS,nMeasTS),dtype=cython_real))
#    cdef np.ndarray[ndim=2,dtype=cython_real_t] Cvv=np.asfortranarray(np.zeros((nact,nact),dtype=cython_real))
#    cdef np.ndarray[ndim=2,dtype=cython_real_t] tmp=np.asfortranarray(np.zeros((nact,nMeasTS),dtype=cython_real))
#
#    compute_Cee_Cvv(nMeas,nMeasTS,nact,<real_t*>Cmm.data,nMeas,<real_t*>cpp.data,nMeasTS,<real_t*>cpm.data,nMeasTS,<real_t*>r.data,nMeasTS,<real_t*>dx.data,nact,<real_t*>Cee.data,nMeasTS,<real_t*>Cvv.data,nact,<real_t*>tmp.data,nact)
#
#    return Cee,Cvv


#
# wrappers for blas/lapack routine
#

#def dtrsm(M,N,np.ndarray[ndim=2,dtype=cython_real_t]A,lda,np.ndarray[ndim=2,dtype=cython_real_t]B,ldb):
#
#    #swap lda,ldb: memory layout
#    cdef real_t *Adata=<real_t *>A.data
#    cdef real_t *Bdata=<real_t *>B.data
#    _dtrsm(M,N,Adata,lda,Bdata,ldb)
#    #A *X = alpha*B
#    #dtrsm(A.shape[0],A.shape[1],A,A.shape[1],B,B.shape[1])
#
def _py_Xtrsm(side,uplo,trans,diag,M,N,alpha,np.ndarray[ndim=2,dtype=cython_real_t]A,np.ndarray[ndim=2,dtype=cython_real_t]B):
    cdef real_t *Adata=<real_t *>A.data
    cdef real_t *Bdata=<real_t *>B.data

    cdef CBLAS_ORDER        layout
    cdef CBLAS_SIDE         s=CblasRight
    cdef CBLAS_UPLO         u=CblasUpper
    cdef CBLAS_TRANSPOSE    t=CblasTrans
    cdef CBLAS_DIAG         d=CblasUnit
    cdef int                lda,ldb

    if(side.upper()[0]=='L'):
        s=CblasLeft
    if(uplo.upper()[0]=='L'):
        u=CblasLower
    if(trans.upper()[0]=='N'):
        t=CblasNoTrans
    if(diag.upper()[0]=='N'):
        d=CblasNonUnit

    if(np.isfortran(A)==np.isfortran(B)):
        if(np.isfortran(A)):
            layout=CblasColMajor
            lda=A.shape[0]
            ldb=B.shape[0]
        else:
            layout=CblasRowMajor
            lda=A.shape[1]
            ldb=B.shape[1]

        try:
            cblas_Xtrsm(layout,s,u,t,d,M,N,alpha,Adata,lda,Bdata,ldb)
        except Exception, e:
            #TODO cannot catch exception: wrong M,lda or N,ldb
            print "Unexpected error:", e


    else:
        print "not the same layout for the arrays"
        return

#
#def dsyr2k(M,N,alpha,np.ndarray[ndim=2,dtype=cython_real_t]A,lda,np.ndarray[ndim=2,dtype=cython_real_t]B,ldb,beta,np.ndarray[ndim=2,dtype=cython_real_t]C,ldc):
#
#    #swap N and M
#    cdef real_t *Adata=<real_t *>A.data
#    cdef real_t *Bdata=<real_t *>B.data
#    cdef real_t *Cdata=<real_t *>C.data
#    _dsyr2k(M,N,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)
#
#    #alpha*A*B' + alpha*B*A' + beta*C
#    #dsyrk(C.shape[0],A.shape[0],alpha,A,A.shape[1],B,B.shape[1],beta,C,C.shape[1])
#
def _py_Xsyr2k(uplo,trans,N,K,alpha,np.ndarray[ndim=2,dtype=cython_real_t] A,np.ndarray[ndim=2,dtype=cython_real_t] B,beta,np.ndarray[ndim=2,dtype=cython_real_t]C):
    cdef real_t *Adata=<real_t*>A.data
    cdef real_t *Bdata=<real_t*>B.data
    cdef real_t *Cdata=<real_t*>C.data

    cdef CBLAS_ORDER        layout
    cdef CBLAS_UPLO         u=CblasUpper
    cdef CBLAS_TRANSPOSE    t=CblasTrans
    cdef int                lda,ldb,ldc

    if(uplo.upper()[0]=='L'):
        u=CblasLower
    if(trans.upper()[0]=='N'):
        t=CblasNoTrans
    if(np.isfortran(A)==np.isfortran(B)):
        if(np.isfortran(A)):
            layout=CblasColMajor
            lda=A.shape[0]
            ldb=B.shape[0]
            ldc=C.shape[0]
        else:
            layout=CblasRowMajor
            lda=A.shape[1]
            ldb=B.shape[1]
            ldc=C.shape[1]

        print layout,uplo,trans,N,K,alpha,'A',lda,'B',ldb,beta,'C',ldc
        #cblas_dsyr2k(layout,u,t,N,K,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)
        cblas_Xsyr2k(CblasRowMajor,CblasUpper,CblasNoTrans,N,K,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)


#def dsymm(M,N,alpha,np.ndarray[ndim=2,dtype=cython_real_t]A,lda,np.ndarray[ndim=2,dtype=cython_real_t]B,ldb,beta,np.ndarray[ndim=2,dtype=cython_real_t]C,ldc):
#
#    cdef real_t *Adata=<real_t *>A.data
#    cdef real_t *Bdata=<real_t *>B.data
#    cdef real_t *Cdata=<real_t *>C.data
#
#    _dsymm(M,N,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)
#    #c := alpha*B*A + beta*C
#    #dsymm(C.shape[0],C.shape[1],alpha,A,A.shape[1],B,B.shape[1],beta,C,C.shape[1]))
#
def _py_Xsymm(side,uplo,M,N,alpha,np.ndarray[ndim=2,dtype=cython_real_t] A, np.ndarray[ndim=2,dtype=cython_real_t]B,beta, np.ndarray[ndim=2,dtype=cython_real_t]C):
    cdef real_t *Adata=<real_t*>A.data
    cdef real_t *Bdata=<real_t*>B.data
    cdef real_t *Cdata=<real_t*>C.data

    cdef CBLAS_ORDER        layout
    cdef CBLAS_SIDE         s=CblasRight
    cdef CBLAS_UPLO         u=CblasUpper
    cdef int                lda,ldb,ldc

    if(side.upper()[0]=='L'):
        s=CblasLeft
    if(uplo.upper()[0]=='L'):
        u=CblasLower

    if(np.isfortran(A)==np.isfortran(B)):
        if(np.isfortran(A)):
            layout=CblasColMajor
            lda=A.shape[0]
            ldb=B.shape[0]
            ldc=C.shape[0]
        else:
            layout=CblasRowMajor
            lda=A.shape[1]
            ldb=B.shape[1]
            ldc=C.shape[1]

        cblas_Xsymm(layout,s,u,M,N,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)

#
#
#def dlacpy(M,N,np.ndarray[ndim=2,dtype=cython_real_t]A,lda,np.ndarray[ndim=2,dtype=cython_real_t]B,ldb):
#
#    cdef real_t *Adata=<real_t *>A.data
#    cdef real_t *Bdata=<real_t *>B.data
#
#    _dlacpy(M,N,Adata,lda,Bdata,ldb)
#    #cpy A in B
#    #dlacpy(A.shape[0],A.shape[1],A,A.shape[1],B,B.shape[1])
#
def _py_Xlacpy(uplo,M,N,np.ndarray[ndim=2,dtype=cython_real_t] A,np.ndarray[ndim=2,dtype=cython_real_t] B):
    cdef real_t *Adata=<real_t*>A.data
    cdef real_t *Bdata=<real_t*>B.data

    #cdef CBLAS_UPLO         u=CblasUpper
    cdef char u='A'
    cdef int                lda,ldb,ldc

    if(uplo.upper()[0]=='U'):
        u='U'
    if(uplo.upper()[0]=='L'):
        u='L'

    if(np.isfortran(A)==np.isfortran(B)):
        if(np.isfortran(A)):
            layout=LAPACK_COL_MAJOR
            lda=A.shape[0]
            ldb=B.shape[0]
        else:
            layout=LAPACK_ROW_MAJOR
            lda=A.shape[1]
            ldb=B.shape[1]


        LAPACKE_Xlacpy(layout,u,M,N,Adata,lda,Bdata,ldb)



#def dgemm(M,N,K,alpha,np.ndarray[ndim=2,dtype=cython_real_t]A,lda,np.ndarray[ndim=2,dtype=cython_real_t]B,ldb,beta,np.ndarray[ndim=2,dtype=cython_real_t]C,ldc):
#    
#    cdef real_t *Adata=<real_t *>A.data
#    cdef real_t *Bdata=<real_t *>B.data
#    cdef real_t *Cdata=<real_t *>C.data
#
#
#    _dgemm(M,N,K,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)
#    #dgemm(C.shape[0],C.shape[1],A.shape[1],2,A,A.shape[1],B,B.shape[1],1,C,C.shape[1])


def _py_Xgemm(transa,transb,M,N,K,alpha,np.ndarray[ndim=2,dtype=cython_real_t] A,np.ndarray[ndim=2,dtype=cython_real_t] B,beta,np.ndarray[ndim=2,dtype=cython_real_t] C):

    cdef real_t *Adata=<real_t*>A.data
    cdef real_t *Bdata=<real_t*>B.data
    cdef real_t *Cdata=<real_t*>C.data

    cdef CBLAS_ORDER        layout
    cdef CBLAS_TRANSPOSE    ta=CblasTrans
    cdef CBLAS_TRANSPOSE    tb=CblasTrans
    cdef int                lda,ldb,ldc


    if(transa.upper()[0]=='N'):
        ta=CblasNoTrans
    if(transb.upper()[0]=='N'):
        tb=CblasNoTrans

    if(np.isfortran(A)==np.isfortran(B)):
        if(np.isfortran(A)):
            layout=CblasColMajor
            lda=A.shape[0]
            ldb=B.shape[0]
            ldc=C.shape[0]
        else:
            layout=CblasRowMajor
            lda=A.shape[1]
            ldb=B.shape[1]
            ldc=C.shape[1]

        cblas_Xgemm(layout,ta,tb,M,N,K,alpha,Adata,lda,Bdata,ldb,beta,Cdata,ldc)



#def py_add(trans,M,N,alpha,np.ndarray[ndim=2,dtype=cython_real_t] A,beta,np.ndarray[ndim=2,dtype=cython_real_t] B):
#    cdef char u
##uplo.upper()[0]
#    cdef char t
##trans.upper()[0]
#
#    if(trans.upper()[0]=='N'):
#        t='N'
#    elif(trans.upper()[0]=='T'):
#        t='T'
#    else:
#        print "wrong trans, expect 'N' or 'T'"
#        return
#
#    Xgeadd(t,M,N,alpha,<real_t*>A.data,A.shape[1],beta,<real_t*>B.data,B.shape[1])


#replace X with s for single precision , d for double
IF USE_SINGLE:
    py_strsm =_py_Xtrsm
    py_ssyr2k=_py_Xsyr2k
    py_ssymm =_py_Xsymm
    py_slacpy=_py_Xlacpy
    py_sgemm =_py_Xgemm
ELSE:
    py_dtrsm =_py_Xtrsm
    py_dsyr2k=_py_Xsyr2k
    py_dsymm =_py_Xsymm
    py_dlacpy=_py_Xlacpy
    py_dgemm =_py_Xgemm
