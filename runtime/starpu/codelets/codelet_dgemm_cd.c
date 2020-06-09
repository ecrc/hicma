/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_zsyrk.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
#include "hicma.h"
#include "hicma_common.h"
#include "auxdescutil.h"
#include "coreblas.h"
#include "coreblas/lapacke.h"
#include <sys/time.h>
#include "runtime/starpu/runtime_codelets.h"
DCODELETS_HEADER(gemmcd_hcore)
#include "hcore_z.h"
#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern int print_index;
extern int print_index_end;
extern int print_mat;
extern void _printmat(double * A, int64_t m, int64_t n, int64_t ld);
extern int use_scratch;
/**
 *
 * @ingroup CORE_HICMA_Complex64_t
 *
 **/
//#if defined(HICMA_DOUBLE)
void HICMA_TASK_dgemm_cd(const MORSE_option_t *options,
                         int n, int nb,
                         double alpha,
                         const MORSE_desc_t *AUV, int ldauv,
                         const MORSE_desc_t *Ark,
                         int Am, int An,
                         const MORSE_desc_t *BUV, int ldbuv,
                         const MORSE_desc_t *Brk,
                         int Bm, int Bn,
                         double beta,
                         const MORSE_desc_t *CD, int ldcd,
                         int Cm, int Cn)
{
    int nAUV = AUV->nb;
    (void)nb;
    struct starpu_codelet *codelet = &cl_dgemmcd_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zsyrk_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);
    
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(AUV, Am, An);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    MORSE_ACCESS_R(Ark, Am, An);
#endif
    MORSE_ACCESS_RW(CD, Cm, Cn);
    MORSE_END_ACCESS_DECLARATION;
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,         &n,                        sizeof(int),
                       STARPU_VALUE,     &alpha,         sizeof(double),
                       STARPU_R,                 RTBLKADDR(AUV, double, Am, An),
                       STARPU_VALUE,       &ldauv,                        sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
                       STARPU_R,                 RTBLKADDR(Ark, double, Am, An),
#endif
                       STARPU_R,                 RTBLKADDR(BUV, double, Bm, Bn),
                       STARPU_VALUE,       &ldbuv,                        sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
                       STARPU_R,                 RTBLKADDR(Brk, double, Bm, Bn),
#endif
                       STARPU_VALUE,      &beta,         sizeof(double),
                       STARPU_RW,                 RTBLKADDR(CD, double, Cm, Cn),
                       STARPU_VALUE,       &ldcd,                        sizeof(int),
                       STARPU_VALUE,    &Am,                 sizeof(int),
                       STARPU_VALUE,    &An,                 sizeof(int),
                       STARPU_VALUE,    &Cm,                 sizeof(int),
                       STARPU_VALUE,    &Cn,                 sizeof(int),
                       STARPU_VALUE,    &nAUV,               sizeof(int),
                       STARPU_SCRATCH,   options->ws_worker,
                       STARPU_VALUE,    &h_work,            sizeof(MORSE_starpu_ws_t *),
                       STARPU_PRIORITY,    options->priority,
                       STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_dgemmcd",
#endif
                       0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_dgemm_cd_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    int n;
    double alpha;
    double *AUV;
    double *Ark;
    double *BUV;
    double *Brk;
    int ldauv;
    int ldbuv;
    double beta;
    double *CD;
    int ldcd;
    int Am, An, Cm, Cn;
    int nAUV;
    
    int idescr = 0;
    AUV = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    if(global_always_fixed_rank == 1){
        fprintf(stderr, "global_always_fixed_rank is one. But HICMA_ALWAYS_FIX_RANK is not defined. Exiting...\n");
        exit(1);
    }
#else
    if(global_always_fixed_rank != 1){
        fprintf(stderr, "global_always_fixed_rank must be one. But it is %d. Exiting...\n", global_always_fixed_rank);
        exit(1);
    }
#endif
    int _Ark;
    if(global_always_fixed_rank == 1){
        _Ark = global_fixed_rank;
    } else {
        _Ark = Ark[0];
    }
    
    BUV = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Brk = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    if(global_always_fixed_rank == 1){
        fprintf(stderr, "global_always_fixed_rank is one. But HICMA_ALWAYS_FIX_RANK is not defined. Exiting...\n");
        exit(1);
    }
#else
    if(global_always_fixed_rank != 1){
        fprintf(stderr, "global_always_fixed_rank must be one. But it is %d. Exiting...\n", global_always_fixed_rank);
        exit(1);
    }
#endif
    int _Brk;
    if(global_always_fixed_rank == 1){
        _Brk = global_fixed_rank;
    } else {
        _Brk = Brk[0];
    }
    CD = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    double* work = NULL;
    work = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    MORSE_starpu_ws_t *h_work;
    starpu_codelet_unpack_args(cl_arg, &n, &alpha,   &ldauv, &ldbuv, &beta, &ldcd, &Am, &An, &Cm, &Cn, &nAUV, &h_work);
    double *AU = AUV;
    int nAU = nAUV/2;
    size_t nelm_AU = (size_t)ldauv * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);
    int ldau = ldauv;
    int ldav = ldauv;
    
    double *BU = BUV;
    int nBU = nAUV/2;
    size_t nelm_BU = (size_t)ldbuv * (size_t)nBU;
    double *BV = &(BUV[nelm_BU]);
    if(HICMA_get_print_index()){
        printf("%d+GEMMCD\t|CD(%d,%d) AUV(%d,%d)%d BUV(%d,%d)%d  N:%d\n",MORSE_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, An, Am, _Brk, n);
    }
    if(print_mat){
        printf("%d\tdgemm_cd-input\n");
        printf("AU:\n");_printmat(AU, n, _Ark, ldau);
        printf("AV:\n");_printmat(AV, ldau, _Ark, ldau);
        printf("BU:\n");_printmat(BU, n, _Brk, ldau);
        printf("BV:\n");_printmat(BV, ldau, _Brk, ldau);
        printf("C:\n");_printmat(CD, n, n, ldcd);
    }
    
    int64_t A_colfactor_ncols = n;
    /// C = C + alpha * A * A'
    /// C = C + alpha * ( (A^u * (A^v * A^v^T) ) * A^u^T)
    /// A^v * B^v^T
    int64_t bufmtx_nrows = _Ark;
    int64_t bufmtx_ncols = _Brk;
    size_t bufmtx_nelm = bufmtx_nrows * bufmtx_ncols;
    
    double* bufmtx = NULL;
    if(use_scratch){
        bufmtx = work;
    } else {
        bufmtx = malloc(bufmtx_nelm * sizeof(double));
    }
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, _Ark, _Brk, n, 1.0, AV, ldau, BU, ldau, 0.0, bufmtx, bufmtx_nrows);
    
    double* bufmtx2 = NULL;
    int64_t bufmtx2_nrows = _Ark;
    int64_t bufmtx2_ncols = n;
    size_t  bufmtx2_nelm = bufmtx2_nrows * bufmtx2_ncols;
    if(use_scratch){
        bufmtx2 = work + bufmtx_nelm;
    } else {
        bufmtx2 = malloc(bufmtx2_nelm * sizeof(double));
    }
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,_Ark,n, _Brk , 1.0, bufmtx, bufmtx_nrows, BV, ldau, 0.0, bufmtx2, bufmtx2_nrows);
    
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, _Ark, -1.0, AU, ldau, bufmtx2, bufmtx2_nrows, 1.0, CD, ldcd);
    
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-GEMMCD\t|CD(%d,%d) AUV(%d,%d)%d N:%d LDA:%d LDCD:%d\t\t\t\t\tSYRK:%.4f\n",MORSE_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, n,
               ldauv, ldcd,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
               );
    }
    if(print_mat){
        printf("%d\tdgemm_cd-output\n");
        _printmat(CD, n, n, ldcd);
    }
    if(use_scratch == 0){
        free(bufmtx);
        free(bufmtx2);
    }
}
#endif /* !defined(MORSE_SIMULATION) */
/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(dgemmcd_hcore, 4, cl_dgemm_cd_hcore_cpu_func)
// CODELETS(zsyrk_hcore, 3, cl_zsyrk_hcore_cpu_func, cl_zsyrk_hcore_cuda_func, STARPU_CUDA_ASYNC)
#else
CODELETS_CPU(dgemmcd_hcore, 6, cl_dgemm_cd_hcore_cpu_func)
// CODELETS(zsyrk_hcore, 4, cl_zsyrk_hcore_cpu_func, cl_zsyrk_hcore_cuda_func, STARPU_CUDA_ASYNC)
#endif
//#endif
