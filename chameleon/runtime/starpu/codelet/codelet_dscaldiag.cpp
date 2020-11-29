/*! @Copyright (c) 2017, King Abdullah University of Science and Technology (KAUST)
 * and Observatoire de Paris Meudon (OBSPM)
 * All rights reserved.
 *
 * MOAO is a software package provided by KAUST and OBSPM
 * @version 0.1.0
 **/

#include "codelet_dscaldiag.hpp"
#include <starpu.h>
#include "chameleon_templates.hpp"
#include <type_traits>

static struct starpu_codelet cl_sscaldiag =
{
    where : STARPU_CPU,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {cl_scaldiag_cpu_func<float>}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {},         //starpu_c<float>uda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    
    nbuffers : 0            //int nbuffers;
};

static struct starpu_codelet cl_dscaldiag =
{
    where : STARPU_CPU,
    can_execute : NULL,     //int (*can_execute)(unsigned workerid, struct starpu_task *task, unsigned nimpl);
    type : STARPU_SEQ,      //enum starpu_codelet_type type;
    max_parallelism : 0,    //int max_parallelism;

    cpu_func : {},          //starpu_cpu_func_t cpu_func STARPU_DEPRECATED;
    cuda_func:{},           //starpu_cuda_func_t cuda_func STARPU_DEPRECATED;
    opencl_func:{},         //starpu_opencl_func_t opencl_func STARPU_DEPRECATED;

    cpu_funcs : {cl_scaldiag_cpu_func<double>}, //starpu_cpu_func_t cpu_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_funcs: {},         //starpu_cuda_func_t cuda_funcs[STARPU_MAXIMPLEMENTATIONS];
    cuda_flags: {},         //char cuda_flags[STARPU_MAXIMPLEMENTATIONS];
    opencl_funcs : {},      //starpu_opencl_func_t opencl_funcs[STARPU_MAXIMPLEMENTATIONS];
    opencl_flags : {},      //char opencl_flags[STARPU_MAXIMPLEMENTATIONS];
    mic_funcs : {},         //starpu_mic_func_t mic_funcs[STARPU_MAXIMPLEMENTATIONS];
    //mpi_ms_funcs :{},       //starpu_mpi_ms_func_t mpi_ms_funcs[STARPU_MAXIMPLEMENTATIONS];
    scc_funcs : {},         //starpu_scc_func_t scc_funcs[STARPU_MAXIMPLEMENTATIONS];
    cpu_funcs_name :{},     //const char *cpu_funcs_name[STARPU_MAXIMPLEMENTATIONS];
    
    nbuffers : 0            //int nbuffers;
};


template<typename T>
int MOAO_STARPU::TASK_Xscaldiag(int M,T alpha,MORSE_desc_t *descA, int k, int lda){

    struct starpu_codelet *cl;
    if(std::is_same<T,float>::value){
        cl=&cl_sscaldiag;
    }
    else{
        cl=&cl_dscaldiag;
    }

    T *data=(T*)descA->get_blkaddr(descA,k,k);
    starpu_insert_task(cl,
            STARPU_VALUE, &alpha,  sizeof(T),
            STARPU_VALUE, &M, sizeof(int),
            STARPU_VALUE, &data,   sizeof(T*),
            STARPU_VALUE, &lda,   sizeof(int),
            0);
}
template int MOAO_STARPU::TASK_Xscaldiag(int M,float alpha,MORSE_desc_t *descA, int k, int lda);
template int MOAO_STARPU::TASK_Xscaldiag(int M,double alpha,MORSE_desc_t *descA, int k, int lda);


template<typename T>
void cl_scaldiag_cpu_func(void *buffers[],void *cl_arg){
    int M, LDA; 
    T *A, alpha;
      
    starpu_codelet_unpack_args(cl_arg, &alpha, &M, &A, &LDA);
    MORSE<T>::cblas_Xscal( M, alpha, A, LDA+1 );
}
template void cl_scaldiag_cpu_func<float>(void *buffers[],void *cl_arg);
template void cl_scaldiag_cpu_func<double>(void *buffers[],void *cl_arg);
