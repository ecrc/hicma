/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dsyrk.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/

#include <sys/time.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(syrk_hcore)

//UPDATE this definition. I only copy-paste from runtime/starpu/codelets/codelet_dcallback.c
/*HICMA_CHAM_CL_CB(zsyrk_hcore,         starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                      ( 1.+   M)*M*N)*/


#include "hcore_d.h"
extern flop_counter counters[FLOP_NUMTHREADS];

int global_always_fixed_rank = 0;
int global_fixed_rank = 0;
int print_index = 0;
int print_index_end = 0;
int print_mat = 0;
extern void _printmat(double * A, int64_t m, int64_t n, int64_t ld);
/**
 *
 * @ingroup CORE_HICMA_Complex64_t
 *
 **/

void HICMA_TASK_dsyrk(const HICMA_option_t *options,
                      HICMA_enum uplo, HICMA_enum trans,
                      int n, int nb,
                      double alpha,
                      const HICMA_desc_t *AUV, int ldauv,
                      const HICMA_desc_t *Ark,
                      int Am, int An,
                      double beta,
                      const HICMA_desc_t *CD, int ldcd,
                      int Cm, int Cn)
{
    int nAUV = AUV->nb;
    (void)nb;
    struct starpu_codelet *codelet = &cl_dsyrk_hcore;
    void (*callback)(void*) =  NULL;
    HICMA_starpu_ws_t *h_work = (HICMA_starpu_ws_t*)(options->ws_host);

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_R(AUV, Am, An);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    HICMA_ACCESS_R(Ark, Am, An);
#endif
    HICMA_ACCESS_RW(CD, Cm, Cn);
    HICMA_END_ACCESS_DECLARATION;
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,      &uplo,                sizeof(HICMA_enum),
            STARPU_VALUE,     &trans,                sizeof(HICMA_enum),
            STARPU_VALUE,         &n,                        sizeof(int),
            STARPU_VALUE,     &alpha,         sizeof(double),
            STARPU_R,                 RTBLKADDR(AUV, double, Am, An),
            STARPU_VALUE,       &ldauv,                        sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
            STARPU_R,                 RTBLKADDR(Ark, double, Am, An),
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
            STARPU_VALUE,    &h_work,            sizeof(HICMA_starpu_ws_t *),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_dsyrk",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_dsyrk_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    HICMA_enum uplo;
    HICMA_enum trans;
    int n;
    double alpha;
    double *AUV;
    double *Ark;
    int ldauv;
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
    CD = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    double* work = NULL;
    work = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    HICMA_starpu_ws_t *h_work;
    starpu_codelet_unpack_args(cl_arg, &uplo, &trans, &n, &alpha,  &ldauv, &beta, &ldcd, &Am, &An, &Cm, &Cn, &nAUV, &h_work);
    double *AU = AUV;
    int nAU = nAUV/2;
    size_t nelm_AU = (size_t)ldauv * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);
    int ldau = ldauv;
    int ldav = ldauv;

    if(print_index){
        printf("%d+SYRK\t|CD(%d,%d) AUV(%d,%d)%d N:%d\n",HICMA_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, n);
    }
    if(print_mat){
        printf("%d\tsyrk-input\n");
        _printmat(AU, n, _Ark, ldau);
        _printmat(AV, ldau, _Ark, ldau);
        _printmat(CD, n, n, ldcd);
    }
    flop_counter flops;
    flops.syrk = 0;
    HCORE_dsyrk(uplo, trans,
        n, _Ark,
        alpha,
        AU, ldau,
        AV, ldav,
        beta,
        CD, ldcd, work,
        &flops
        );
    int myid = HICMA_RUNTIME_thread_rank(NULL);
    counters[myid].syrk += flops.syrk; 

    /*cblas_dsyrk(*/
        /*CblasColMajor,*/
        /*(CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,*/
        /*n, k,*/
        /*CBLAS_SADDR(alpha), A, lda,*/
        /*CBLAS_SADDR(beta), C, ldc);*/
    if(print_index || print_index_end){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-SYRK\t|CD(%d,%d) AUV(%d,%d)%d N:%d LDA:%d LDCD:%d\t\t\t\t\tSYRK:%.4f\n",HICMA_My_Mpi_Rank(),Cm, Cn, Am, An, _Ark, n,
                ldauv, ldcd,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
    if(print_mat){
        printf("%d\tsyrk-output\n");
        _printmat(CD, n, n, ldcd);
    }
}
#endif /* !defined(HICMA_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(dsyrk_hcore, 3, cl_dsyrk_hcore_cpu_func)
#else
CODELETS_CPU(dsyrk_hcore, 4, cl_dsyrk_hcore_cpu_func)
#endif
