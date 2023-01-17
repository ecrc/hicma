/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zgemm_cd.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/

#include <sys/time.h>
#include <runtime/starpu/hicma_starpu.h>
#include <hicma.h>
#include <misc/auxdescutil.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

ZCODELETS_HEADER(gemmcd_hcore)
#include "hcore_z.h"
#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

int isPercentageused =0;
extern float percent1;
float percent2 =0.5;
extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern int print_index;
extern int print_index_end;
extern int print_mat;
extern void _printmat_complex_complex(HICMA_Complex64_t * A, int64_t m, int64_t n, int64_t ld);
extern int use_scratch;

/**
 *
 * @ingroup CORE_HICMA_Complex64_t
 *
 **/

//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zgemm_cd(const HICMA_option_t *options,
                         int n, int nb,
                         HICMA_Complex64_t alpha,
                         const HICMA_desc_t *AUV, int ldauv,
                         const HICMA_desc_t *Ark,
                         int Am, int An,
                         const HICMA_desc_t *BUV, int ldbuv,
                         const HICMA_desc_t *Brk,
                         int Bm, int Bn,
                         HICMA_Complex64_t beta,
                         const HICMA_desc_t *CD, int ldcd,
                         int Cm, int Cn)
{
    int nAUV = AUV->nb;
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgemmcd_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zsyrk_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    HICMA_starpu_ws_t *h_work = (HICMA_starpu_ws_t*)(options->ws_host);
    
    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_R(AUV, Am, An);
    HICMA_ACCESS_R(BUV, Bm, Bn);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    HICMA_ACCESS_R(Ark, Am, An);
    HICMA_ACCESS_R(Brk, Bm, Bn);
#endif
    //printf("%s %d: A:%d,%d B:%d,%d C:%d,%d ldcd:%d\n", __FILE__, __LINE__, Am, An, Bm, Bn, Cm, Cn, ldcd); fflush(stdout);
    HICMA_ACCESS_RW(CD, Cm, Cn);
    HICMA_END_ACCESS_DECLARATION;
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,         &n,                        sizeof(int),
                       STARPU_VALUE,     &alpha,         sizeof(HICMA_Complex64_t),
                       STARPU_R,                 RTBLKADDR(AUV, HICMA_Complex64_t, Am, An),
                       STARPU_VALUE,       &ldauv,                        sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
                       STARPU_R,                 RTBLKADDR(Ark, double, Am, An),
#endif
                       STARPU_R,                 RTBLKADDR(BUV, HICMA_Complex64_t, Bm, Bn),
                       STARPU_VALUE,       &ldbuv,                        sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
                       STARPU_R,                 RTBLKADDR(Brk, double, Bm, Bn),
#endif
                       STARPU_VALUE,      &beta,         sizeof(HICMA_Complex64_t),
                       STARPU_RW,                 RTBLKADDR(CD, HICMA_Complex64_t, Cm, Cn),
                       STARPU_VALUE,       &ldcd,                        sizeof(int),
                       STARPU_VALUE,    &Am,                 sizeof(int),
                       STARPU_VALUE,    &An,                 sizeof(int),
                       STARPU_VALUE,    &Cm,                 sizeof(int),
                       STARPU_VALUE,    &Cn,                 sizeof(int),
                       STARPU_VALUE,    &nAUV,               sizeof(int),
                       STARPU_SCRATCH,   options->ws_worker,
                       STARPU_VALUE,    &h_work,            sizeof(HICMA_starpu_ws_t *),
                       STARPU_PRIORITY,    options->priority,
                       STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zgemmcd",
#endif
                       0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgemm_cd_hcore_cpu_func(void *descr[], void *cl_arg)
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
    HICMA_Complex64_t alpha;
    HICMA_Complex64_t *AUV;
    double *Ark;
    HICMA_Complex64_t *BUV;
    double *Brk;
    int ldauv;
    int ldbuv;
    HICMA_Complex64_t beta;
    HICMA_Complex64_t *CD;
    int ldcd;
    int Am, An, Cm, Cn;
    int nAUV;
    
    int idescr = 0;
    AUV = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Ark = (double*)STARPU_MATRIX_GET_PTR(descr[idescr++]);
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
    
    BUV = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
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
    CD = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    HICMA_Complex64_t* work = NULL;
    work = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    HICMA_starpu_ws_t *h_work;
    starpu_codelet_unpack_args(cl_arg, &n, &alpha,   &ldauv, &ldbuv, &beta, &ldcd, &Am, &An, &Cm, &Cn, &nAUV, &h_work);

    /*split AU and AV*/
    HICMA_Complex64_t *AU = AUV;
   int nAU, nBU;

    nAU = nAUV/2;
    size_t nelm_AU = (size_t)ldauv * (size_t)nAU;
    HICMA_Complex64_t *AV = &(AUV[nelm_AU]);
    int ldau = ldauv;
    int ldav = ldauv;
    
   /*split BU and BV */
    HICMA_Complex64_t *BU = BUV;
   if(isPercentageused){
    int spb=percent2*_Brk;
    if (spb==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); spb=1;}
     nBU = _Brk+spb;
    }
    else{ nBU = nAUV/2;}
    size_t nelm_BU = (size_t)ldbuv * (size_t)nBU;
    HICMA_Complex64_t *BV = &(BUV[nelm_BU]);
    if(HICMA_get_print_index()){
        printf("%d+GEMMCD\t|CD(%d,%d) AUV(%d,%d)%d BUV(%d,%d)%d  N:%d\n",HICMA_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, An, Am, _Brk, n);
    }
    if(print_mat){
        printf("%d\tzgemm_cd-input\n");
        printf("AU:\n");_printmat_complex(AU, n, _Ark, ldau);
        printf("AV:\n");_printmat_complex(AV, ldau, _Ark, ldau);
        printf("BU:\n");_printmat_complex(BU, n, _Brk, ldau);
        printf("BV:\n");_printmat_complex(BV, ldau, _Brk, ldau);
        printf("C:\n");_printmat_complex(CD, n, n, ldcd);
    }
    
    int64_t A_colfactor_ncols = n;
    /// C = C + alpha * A * A'
    /// C = C + alpha * ( (A^u * (A^v * A^v^T) ) * A^u^T)
    /// A^v * B^v^T
    int64_t bufmtx_nrows = _Ark;
    int64_t bufmtx_ncols = _Brk;
    size_t bufmtx_nelm = bufmtx_nrows * bufmtx_ncols;
    
    HICMA_Complex64_t* bufmtx = NULL;
    if(use_scratch){
        bufmtx = work;
    } else {
        bufmtx = malloc(bufmtx_nelm * sizeof(HICMA_Complex64_t));
    }
    HICMA_Complex64_t alphacom=1.0, betacom=0.0;
    cblas_zgemm(CblasColMajor, CblasTrans, CblasNoTrans, _Ark, _Brk, n, &alphacom, AV, ldau, BU, ldau, &betacom, bufmtx, bufmtx_nrows);

    HICMA_Complex64_t* bufmtx2 = NULL;
    int64_t bufmtx2_nrows = _Ark;
    int64_t bufmtx2_ncols = n;
    size_t  bufmtx2_nelm = bufmtx2_nrows * bufmtx2_ncols;
    if(use_scratch){
        bufmtx2 = work + bufmtx_nelm;
    } else {
        bufmtx2 = malloc(bufmtx2_nelm * sizeof(HICMA_Complex64_t));
    }
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans,_Ark,n, _Brk , &alphacom, bufmtx, bufmtx_nrows, BV, ldau, &betacom, bufmtx2, bufmtx2_nrows);
    
    alphacom=-1.0; betacom=1.0;
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, _Ark, &alphacom, AU, ldau, bufmtx2, bufmtx2_nrows, &betacom, CD, ldcd);
    
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-GEMMCD\t|CD(%d,%d) AUV(%d,%d)%d N:%d LDA:%d LDCD:%d\t\t\t\t\tSYRK:%.4f\n",HICMA_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, n,
               ldauv, ldcd,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
               );
    }
   // if(print_mat){
        //printf("%d\tzgemm_cd-output\n");
       //_printmat_complex_complex(CD, n, n, ldcd);
   // }
    if(use_scratch == 0){
        free(bufmtx);
        free(bufmtx2);
    }
}
#endif /* !defined(HICMA_SIMULATION) */
/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(zgemmcd_hcore, 4, cl_zgemm_cd_hcore_cpu_func)
// CODELETS(zsyrk_hcore, 3, cl_zsyrk_hcore_cpu_func, cl_zsyrk_hcore_cuda_func, STARPU_CUDA_ASYNC)
#else
CODELETS_CPU(zgemmcd_hcore, 6, cl_zgemm_cd_hcore_cpu_func)
// CODELETS(zsyrk_hcore, 4, cl_zsyrk_hcore_cpu_func, cl_zsyrk_hcore_cuda_func, STARPU_CUDA_ASYNC)
#endif
//#endif
