/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_ztrsm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
#include "hicma.h"
#include "hicma_common.h"
#include "auxdescutil.h"
#include "coreblas.h"
#include "coreblas/lapacke.h"
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
//#include "runtime/starpu/include/runtime_codelet_z.h"

#include <sys/time.h>

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(trsm_hcore)

//UPDATE this definition. I only copy-paste from runtime/starpu/codelets/codelet_zcallback.c
/*CHAMELEON_CL_CB(ztrsm_hcore,         starpu_matrix_get_nx(task->handles[1]), starpu_matrix_get_ny(task->handles[1]), 0,                                               M*M*N)*/

#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

int trsm_print_index_end = 0;

void HICMA_TASK_ztrsm(const MORSE_option_t *options,
                      MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
                      int m, 
                      double alpha, const MORSE_desc_t *A, int Am, int An, int lda,
                      const MORSE_desc_t *BUV, int Bm, int Bn, int ldb, const MORSE_desc_t *Brk)
{
    int nBUV = BUV->nb;
    struct starpu_codelet *codelet = &cl_ztrsm_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_ztrsm_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    int sizeA = lda*m;
    int sizeB = ldb; //*nb; //@KADIR converted n to nb FIXME Size of B will be determined at runtime!!!
    int execution_rank = BUV->get_rankof( BUV, Bm, Bn );
    int rank_changed=0;
    (void)execution_rank;

    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char* env = getenv("MORSE_COMM_FACTOR_THRESHOLD");
    if (env != NULL)
        threshold = (unsigned)atoi(env);
    else
        threshold = 10;
    if ( sizeA > threshold*sizeB ){
        execution_rank = A->get_rankof( A, Am, An );
        rank_changed=1;
    }
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(A, Am, An);
    MORSE_ACCESS_RW(BUV, Bm, Bn);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    MORSE_ACCESS_R(Brk, Bm, Bn);
#endif
    if (rank_changed)
        MORSE_RANK_CHANGED(execution_rank);
    MORSE_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &side,               sizeof(MORSE_enum),
            STARPU_VALUE,    &uplo,               sizeof(MORSE_enum),
            STARPU_VALUE,    &transA,             sizeof(MORSE_enum),
            STARPU_VALUE,    &diag,               sizeof(MORSE_enum),
            STARPU_VALUE,    &m,                  sizeof(int),
            STARPU_VALUE,    &alpha,              sizeof(double),
            STARPU_R,         RTBLKADDR(A, double, Am, An),
            STARPU_VALUE,    &lda,                sizeof(int),
            STARPU_RW,        RTBLKADDR(BUV, double, Bm, Bn),
            STARPU_VALUE,    &ldb,                sizeof(int),
#if !defined(HICMA_ALWAYS_FIX_RANK)
            STARPU_R,        RTBLKADDR(Brk, double, Bm, Bn),
#endif
            STARPU_VALUE,    &Am,                 sizeof(int),
            STARPU_VALUE,    &An,                 sizeof(int),
            STARPU_VALUE,    &Bm,                 sizeof(int),
            STARPU_VALUE,    &Bn,                 sizeof(int),
            STARPU_VALUE,    &nBUV,               sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
            STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_ztrsm",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrsm_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    MORSE_enum side;
    MORSE_enum uplo;
    MORSE_enum transA;
    MORSE_enum diag;
    int m;
    double alpha;
    double *A;
    int lda;
    double *BUV;
    int ldb;
    double *Brk;
    int Am;
    int An;
    int Bm;
    int Bn;
    int nBUV;

    A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    BUV = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Brk = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    if(HICMA_get_always_fixed_rank() == 1){
        fprintf(stderr, "global_always_fixed_rank is one. But HICMA_ALWAYS_FIX_RANK is not defined. Exiting...\n");
        exit(1);
    }
#else
    if(HICMA_get_always_fixed_rank() != 1){
        fprintf(stderr, "global_always_fixed_rank must be one. But it is %d. Exiting...\n", HICMA_get_always_fixed_rank());
        exit(1);
    }
#endif
    int _Brk;
    if(HICMA_get_always_fixed_rank() == 1){
        _Brk = HICMA_get_fixed_rank();
    } else {
        _Brk = Brk[0];
    }

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &m,  &alpha, &lda, &ldb, &Am, &An, &Bm, &Bn, &nBUV);

    int nBU = nBUV/2;
    size_t nelm_BU = (size_t)ldb * (size_t)nBU;
    double *B = &(BUV[nelm_BU]);

    /*CORE_ztrsm(side, uplo,*/
        /*transA, diag,*/
        /*m, n,*/
        /*alpha, A, lda,*/
        /*B, ldb);*/
    if(HICMA_get_print_index() == 1){
        printf("%d+TRSM\t|AD(%d,%d) BV(%d,%d)%d m:%d lda(11):%d ldb(12):%d\n",MORSE_My_Mpi_Rank(),Am,An, Bm, Bn, _Brk, m, lda, ldb);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\ttrsm-input A\n", __LINE__);
        _printmat(A, m, m, lda);
        printf("%d\ttrsm-input B\n", __LINE__);
        _printmat(B, m, _Brk, ldb);
    }
    cblas_dtrsm(
        CblasColMajor,
        (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
        (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
        m,
        _Brk,
        CBLAS_SADDR(alpha), A, lda,
        B, ldb);
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1 || trsm_print_index_end){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-TRSM\t|AD(%d,%d)%dx%d-%d BV(%d,%d)%dx%d-%d m:%d\t\t\t\tTRSM: %.4f\n",MORSE_My_Mpi_Rank(),Am,An, m, m, lda,Bm, Bn, m, _Brk, ldb, m,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\ttrsm-output\n", __LINE__);
        _printmat(B, m, _Brk, ldb);
    }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(ztrsm_hcore, 2, cl_ztrsm_hcore_cpu_func)
// CODELETS(ztrsm_hcore, 2, cl_ztrsm_hcore_cpu_func, cl_ztrsm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#else
CODELETS_CPU(ztrsm_hcore, 3, cl_ztrsm_hcore_cpu_func)
// CODELETS(ztrsm_hcore, 3, cl_ztrsm_hcore_cpu_func, cl_ztrsm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#endif
