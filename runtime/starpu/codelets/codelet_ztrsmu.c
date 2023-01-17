/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_ztrsmu.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Ali Charara
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

#include <sys/time.h>
#include <hicma.h>
#include <misc/auxdescutil.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

ZCODELETS_HEADER(trsmu_hcore)

//UPDATE this definition. I only copy-paste from runtime/starpu/codelets/codelet_zcallback.c
/*HICMA_CHAM_CL_CB(ztrsm_hcore,         starpu_matrix_get_nx(task->handles[1]), starpu_matrix_get_ny(task->handles[1]), 0,                                               M*M*N)*/


int trsmu_print_index_end = 0;
//#if defined(HICMA_COMPLEX)
#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) &(_val)
void HICMA_TASK_hcore_ztrsmu(const HICMA_option_t *options,
                       HICMA_enum side, HICMA_enum uplo, HICMA_enum transA, HICMA_enum diag,
                       int m,
                       HICMA_Complex64_t alpha, const HICMA_desc_t *A, int Am, int An, int lda,
                       const HICMA_desc_t *BUV, int Bm, int Bn, int ldb, const HICMA_desc_t *Brk)
{
    int nBUV = BUV->nb;
    struct starpu_codelet *codelet = &cl_ztrsmu_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_ztrsm_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    int sizeA = lda*m;
    int sizeB = ldb; //*nb; //@KADIR converted n to nb FIXME Size of B will be determined at runtime!!!
    int execution_rank = BUV->get_rankof( BUV, Bm, Bn );
    int rank_changed=0;
    (void)execution_rank;
    
    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char* env = getenv("HiCMA_COMM_FACTOR_THRESHOLD");
    if (env != NULL)
        threshold = (unsigned)atoi(env);
    else
        threshold = 10;
    if ( sizeA > threshold*sizeB ){
        execution_rank = A->get_rankof( A, Am, An );
        rank_changed=1;
    }
    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_R(A, Am, An);
    HICMA_ACCESS_RW(BUV, Bm, Bn);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    HICMA_ACCESS_R(Brk, Bm, Bn);
#endif
    if (rank_changed)
        HICMA_RANK_CHANGED(execution_rank);
    HICMA_END_ACCESS_DECLARATION;
    
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,    &side,               sizeof(HICMA_enum),
                       STARPU_VALUE,    &uplo,               sizeof(HICMA_enum),
                       STARPU_VALUE,    &transA,             sizeof(HICMA_enum),
                       STARPU_VALUE,    &diag,               sizeof(HICMA_enum),
                       STARPU_VALUE,    &m,                  sizeof(int),
                       STARPU_VALUE,    &alpha,              sizeof(HICMA_Complex64_t),
                       STARPU_R,         RTBLKADDR(A, HICMA_Complex64_t, Am, An),
                       STARPU_VALUE,    &lda,                sizeof(int),
                       STARPU_RW,        RTBLKADDR(BUV, HICMA_Complex64_t, Bm, Bn),
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
#if defined(HICMA_USE_MPI)
                       STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_ztrsmu",
#endif
                       0);
}


#if !defined(CHAMELEON_SIMULATION)
static void cl_ztrsmu_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    HICMA_enum side;
    HICMA_enum uplo;
    HICMA_enum transA;
    HICMA_enum diag;
    int m;
    HICMA_Complex64_t alpha;
    HICMA_Complex64_t *A;
    int lda;
    HICMA_Complex64_t *BUV;
    int ldb;
    double *Brk;
    int Am;
    int An;
    int Bm;
    int Bn;
    int nBUV;
    
    A = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    BUV = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Brk = (double*)STARPU_MATRIX_GET_PTR(descr[2]);
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
   int nBU;    
      if(0)printf("%s %s %d: _Brk:%d\n", __FILE__, __func__, __LINE__, _Brk);

    nBU = nBUV/2;
    size_t nelm_BU = (size_t)ldb * (size_t)nBU;
    HICMA_Complex64_t *B = &(BUV[0]);
    
    /*CORE_ztrsm(side, uplo,*/
    /*transA, diag,*/
    /*m, n,*/
    /*alpha, A, lda,*/
    /*B, ldb);*/
    if(HICMA_get_print_index() == 1){
        printf("%d+TRSMU\t|AD(%d,%d) BU(%d,%d)%d m:%d lda(11):%d ldb(12):%d\n",HICMA_My_Mpi_Rank(),Am,An, Bm, Bn, _Brk, m, lda, ldb);
    }
    //if(1 || HICMA_get_print_mat() == 1){
        //printf("%d\ttrsmu-input A\n", __LINE__);
        //_printmat_complex(A, m, m, lda);
       //printf("%d\ttrsmu-input B\n", __LINE__);
       //_printmat_complex(B, m, _Brk, ldb);
    // }
    cblas_ztrsm(
                CblasColMajor,
                (CBLAS_SIDE)side, (CBLAS_UPLO)uplo,
                (CBLAS_TRANSPOSE)transA, (CBLAS_DIAG)diag,
                m,
                _Brk,
                CBLAS_SADDR(alpha), A, lda,
                B, ldb);
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1 || trsmu_print_index_end){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-TRSMU\t|AD(%d,%d)%dx%d-%d BV(%d,%d)%dx%d-%d m:%d\t\t\t\tTRSMu: %.4f\n",HICMA_My_Mpi_Rank(),Am,An, m, m, lda,Bm, Bn, m, _Brk, ldb, m,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
               );
    }
   // if(HICMA_get_print_mat() == 1){
        //printf("%d\ttrsmu-output\n", __LINE__);
        //_printmat_complex(B, m, _Brk, ldb);
  //  }
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(ztrsmu_hcore, 2, cl_ztrsmu_hcore_cpu_func)
// CODELETS(ztrsm_hcore, 2, cl_ztrsm_hcore_cpu_func, cl_ztrsm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#else
CODELETS_CPU(ztrsmu_hcore, 3, cl_ztrsmu_hcore_cpu_func)
// CODELETS(ztrsm_hcore, 3, cl_ztrsm_hcore_cpu_func, cl_ztrsm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#endif
//#endif
