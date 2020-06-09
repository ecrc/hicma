/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zgemm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
#include "morse.h"
#include "hicma.h"
#include "hicma_common.h"
#include "runtime/starpu/chameleon_starpu.h"
//#include "runtime/starpu/include/runtime_codelet_z.h"

#include <sys/time.h>

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(gemm_hcore)


//UPDATE this definition. I only copy-paste from runtime/starpu/codelets/codelet_zcallback.c
/*CHAMELEON_CL_CB(zgemm_hcore,         starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_ny(task->handles[2]), starpu_matrix_get_ny(task->handles[0]),     2. *M*N*K) [> If A^t, computation is wrong <]*/

#include "hcore_z.h"

int codelet_zgemm_print_index = 0;
extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern int print_index;
extern int print_index_end;
extern int print_mat;
extern void _printmat_complex(MORSE_Complex64_t * A, int64_t m, int64_t n, int64_t ld);
/**
 *
 * @ingroup hcore_zgemm
 *
 **/

//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zgemm(const MORSE_option_t *options,
                      MORSE_enum transA, int transB,
                      int m, int n,
                      MORSE_Complex64_t alpha,
                      const MORSE_desc_t *AUV,
                      const MORSE_desc_t *Ark,
                      int Am, int An, int lda,
                      const MORSE_desc_t *BUV,
                      const MORSE_desc_t *Brk,
                      int Bm, int Bn, int ldb,
                      MORSE_Complex64_t beta,
                      const MORSE_desc_t *CUV,
                      const MORSE_desc_t *Crk,
                      int Cm, int Cn, int ldc,
                      int rk,
                      int maxrk,
                      double acc
                      )
{
    int nAUV = AUV->nb;
    int nBUV = BUV->nb;
    int nCUV = CUV->nb;
    struct starpu_codelet *codelet = &cl_zgemm_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zgemm_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);
    /*printf("%s %d:\t%p %p\n", __FILE__, __LINE__, h_work, options->ws_host);*/
    
    int sizeA = lda*nAUV; //FIXME Think about scheduling of tasks according to sizes of the matrices
    int sizeB = ldb*nBUV;
    int sizeC = ldc*nCUV;
    int execution_rank = CUV->get_rankof( CUV, Cm, Cn );
    int rank_changed=0;
    (void)execution_rank;
    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char* env = getenv("MORSE_COMM_FACTOR_THRESHOLD");
    
    int ifval = 0, elseifval = 0, initialval = execution_rank;
    if (env != NULL)
        threshold = (unsigned)atoi(env);
    else
        threshold = 10;
    if ( sizeA > threshold*sizeC ){
        execution_rank = AUV->get_rankof( AUV, Am, An );
        ifval = execution_rank;
        rank_changed = 1;
    }else if( sizeB > threshold*sizeC ){
        execution_rank = BUV->get_rankof( BUV, Bm, Bn );
        elseifval = execution_rank;
        rank_changed = 1;
    }
    //printf("%d,%d %d,%d %d,%d\n", Am, An, Bm, Bn, Cm, Cn);
    //printf("initialval:\t%d if:%d\t else:\t%d rc:\t%d\n", initialval, ifval, elseifval, rank_changed);
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(AUV, Am, An);
    MORSE_ACCESS_R(BUV, Bm, Bn);
    MORSE_ACCESS_RW(CUV, Cm, Cn);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    MORSE_ACCESS_R(Ark, Am, An);
    MORSE_ACCESS_R(Brk, Bm, Bn);
    MORSE_ACCESS_RW(Crk, Cm, Cn);
#endif
    if (rank_changed)
        MORSE_RANK_CHANGED(execution_rank);
    MORSE_END_ACCESS_DECLARATION;
    //printf("%s %d n:%d\n", __func__, __LINE__,n );
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,     &transA,            sizeof(MORSE_enum),
                       STARPU_VALUE,     &transB,            sizeof(MORSE_enum),
                       STARPU_VALUE,     &m,                 sizeof(int),
                       STARPU_VALUE,     &n,                 sizeof(int),
                       STARPU_VALUE,     &alpha,             sizeof(MORSE_Complex64_t),
                       STARPU_R,         RTBLKADDR(AUV, MORSE_Complex64_t, Am, An),
                       STARPU_VALUE,     &lda,               sizeof(int),
                       STARPU_R,         RTBLKADDR(BUV, MORSE_Complex64_t, Bm, Bn),
                       STARPU_VALUE,     &ldb,               sizeof(int),
                       STARPU_VALUE,     &beta,              sizeof(double),
                       STARPU_RW,        RTBLKADDR(CUV, MORSE_Complex64_t, Cm, Cn),
#if !defined(HICMA_ALWAYS_FIX_RANK)
                       STARPU_R,         RTBLKADDR(Ark, double, Am, An),
                       STARPU_R,         RTBLKADDR(Brk, double, Bm, Bn),
                       STARPU_RW,        RTBLKADDR(Crk, double, Cm, Cn),
#endif
                       STARPU_VALUE,     &ldc,               sizeof(int),
                       STARPU_VALUE,     &rk,                 sizeof(int),
                       STARPU_VALUE,     &maxrk,                 sizeof(int),
                       STARPU_VALUE,     &acc,              sizeof(double),
                       STARPU_VALUE,     &Am,                 sizeof(int),
                       STARPU_VALUE,     &An,                 sizeof(int),
                       STARPU_VALUE,     &Bm,                 sizeof(int),
                       STARPU_VALUE,     &Bn,                 sizeof(int),
                       STARPU_VALUE,     &Cm,                 sizeof(int),
                       STARPU_VALUE,     &Cn,                 sizeof(int),
                       STARPU_VALUE,     &nAUV,                 sizeof(int),
                       STARPU_VALUE,     &nBUV,                 sizeof(int),
                       STARPU_VALUE,     &nCUV,                 sizeof(int),
                       STARPU_SCRATCH,   options->ws_worker,
                       STARPU_VALUE,     &h_work,            sizeof(MORSE_starpu_ws_t *),
                       STARPU_PRIORITY,  options->priority,
                       STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
                       STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zgemm",
#endif
                       0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgemm_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    MORSE_Complex64_t alpha;
    MORSE_Complex64_t *AUV = NULL;
    MORSE_Complex64_t *AD = NULL;
    double *Ark = NULL;
    int lda;
    MORSE_Complex64_t *BUV = NULL;
    MORSE_Complex64_t *BD = NULL;
    double *Brk = NULL;
    int ldb;
    MORSE_Complex64_t beta;
    MORSE_Complex64_t *CUV = NULL;
    MORSE_Complex64_t *CD = NULL;
    double *Crk = NULL;
    int ldc;
    int rk;
    int maxrk;
    double acc ;
    int nAUV;
    int nBUV;
    int nCUV;
    
    
    int idescr = 0;
    AUV  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    BUV  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    CUV  = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    Brk = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    Crk = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#else
    double _gemm_rank = global_fixed_rank;
    Ark = &_gemm_rank;
    Brk = &_gemm_rank;
    Crk = &_gemm_rank;
#endif
   
    int _Ark;
    if(global_always_fixed_rank == 1){
        _Ark = global_fixed_rank;
    } else {
        _Ark = Ark[0];
    }

    int _Brk;
    if(global_always_fixed_rank == 1){
        _Brk = global_fixed_rank;
    } else {
        _Brk = Brk[0];
    }

    int  _Crk;
    if(global_always_fixed_rank == 1){
        _Crk = global_fixed_rank;
    } else {
        _Crk = Crk[0];
    } 
    MORSE_Complex64_t* work = NULL;
    work = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    
    int Am, An, Bm, Bn, Cm, Cn;
    
    MORSE_starpu_ws_t *h_work;
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &alpha, &lda, &ldb, &beta, &ldc, &rk, &maxrk, &acc, &Am, &An, &Bm, &Bn, &Cm, &Cn, &nAUV, &nBUV, &nCUV, &h_work);
    
    //printf("%d,%d:%g\t%d,%d:%g\t%d,%d:%g\n", Am, An, *Ark, Bm, Bn, *Brk, Cm, Cn, *Crk);
    MORSE_Complex64_t *AU = AUV;
    MORSE_Complex64_t *BU = BUV;
    MORSE_Complex64_t *CU = CUV;
   int nAU, nBU, nCU;

   if(isPrecentageused){
    int spa=precent2*_Ark;
    if (spa==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); spa=1;}
     nAU = _Ark+spa;
    }
    else { nAU = nAUV/2;}
    if(0)printf("%s %s %d: Am:%d, An:%d _Ark:%d\n", __FILE__, __func__, __LINE__, Am, An, _Ark); 
    size_t nelm_AU = (size_t)lda * (size_t)nAU;
    MORSE_Complex64_t *AV = &(AUV[nelm_AU]);

   if(isPrecentageused){
    int spb=precent2*_Brk;
    if (spb==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); spb=1;}
     nBU = _Brk+spb;
    }
    else{ nBU = nBUV/2;}
     if(0)printf("%s %s %d: Bm:%d, Bn:%d _Brk:%d\n", __FILE__, __func__, __LINE__, Bm, Bn, _Brk);
    size_t nelm_BU = (size_t)ldb * (size_t)nBU;
    MORSE_Complex64_t *BV = &(BUV[nelm_BU]);
   
   if(isPrecentageused){ 
    int spc=precent2*_Crk;
    if (spc==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); spc=1;}
     nCU = _Crk+spc; 
   }
   else{ nCU = nCUV/2;}
    size_t nelm_CU = (size_t)ldc * (size_t)nCU;
    MORSE_Complex64_t *CV = &(CUV[nelm_CU]);
    
    double old_Crk = Crk[0];
    char datebuf_start[128];
    datebuf_start[0] = '\0';
    if( HICMA_get_print_index() || codelet_zgemm_print_index == 1){
        time_t timer;
        struct tm* tm_info;
        gettimeofday (&tvalAfter, NULL);
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf_start, 26, "%Y-%m-%d %H:%M:%S",tm_info); \
        printf("%d+GEMM\t|CUV(%d,%d)%g AUV(%d,%d)%g BUV(%d,%d)%g\t\t\t\t\tGEMM: %s\n",MORSE_My_Mpi_Rank(), Cm, Cn, old_Crk, Am, An, Ark[0], Bm, Bn, Brk[0], datebuf_start);
   }
    
   // if(print_mat){
        //printf("%d\tzgemm-input\n");
        //int _Ark = (int) *Ark;
        //printf("AU:\n");_printmat_complex(AU, n, _Ark, lda);
       // printf("AV:\n");_printmat_complex(AV, n, _Ark, lda);
        //int _Brk = (int) *Brk;
        //printf("BU:\n");_printmat_complex(BU, n, _Brk, ldb);
        //printf("BV:\n");_printmat_complex(BV, n, _Brk, ldb);
    //}
   
     _Crk = (int) *Crk;
    //    printf("OUTCU:\n");_printmat_complex(CU, n, _Crk, ldc);
    //    printf("OUTCV:\n");_printmat_complex(CV, n, _Crk, ldc);
         if(0)printf("%s %s %d: 1_Crk:%d, Cm:%d, Cn:%d\n", __FILE__, __func__, __LINE__, _Crk, Cm, Cn); 
    int isTransA = transA == MorseTrans;
    int isTransB = transB == MorseTrans;
    
    if(HICMA_get_use_fast_hcore_zgemm() == 1){
        HCORE_zgemm_fast(transA, transB,
                         m, n,
                         alpha, (isTransA ? AV : AU), (isTransA ? AU : AV), Ark, lda,
                         (isTransB ? BU : BV), (isTransB ? BV : BU), Brk, ldb,
                         beta, CU, CV, Crk, ldc, rk, maxrk, acc, work);
    }else{
        HCORE_zgemm(transA, transB,
                    m, n,
                    alpha, (isTransA ? AV : AU), (isTransA ? AU : AV), Ark, lda,
                    (isTransB ? BU : BV), (isTransB ? BV : BU), Brk, ldb,
                    beta, CU, CV, Crk, ldc, rk, maxrk, acc, work);
    }
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1 || codelet_zgemm_print_index == 1){
        char datebuf[128];
        time_t timer;
        struct tm* tm_info;
        gettimeofday (&tvalAfter, NULL);
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S",tm_info); \
        printf("%d-GEMM\t|CUV(%d,%d)%g->%g AUV(%d,%d)%g BUV(%d,%d)%g acc:%e rk:%d maxrk:%d\t\t\tGEMM: %.4f\t%s---%s\n",MORSE_My_Mpi_Rank(),Cm, Cn, old_Crk, Crk[0],  Am, An, Ark[0], Bm, Bn, Brk[0],
               acc, rk, maxrk,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0,
               datebuf_start, datebuf
               );
    }

        int _Crk_new = (int) *Crk;
         if(0)printf("%s %s %d: 2_Ark:%d Am:%d, An:%d\n", __FILE__, __func__, __LINE__, _Ark, Am, An);
         if(0)printf("%s %s %d: 2_Brk:%d Bm:%d, Bn:%d\n", __FILE__, __func__, __LINE__, _Brk, Bm, Bn);
         if(0)printf("%s %s %d: 2_Crk:%d Cm:%d, Cn:%d\n", __FILE__, __func__, __LINE__, _Crk_new, Cm, Cn);
        //int _Crk = (int) *Crk;
        //printf("CU:\n");_printmat_complex(CU, n, _Crk_new, ldc);
        //printf("CV:\n");_printmat_complex(CV, n, _Crk_new, ldc);

}
#endif /* !defined(MORSE_SIMULATION) */

/*
 * Codelet definition
 */
#if defined(HICMA_ALWAYS_FIX_RANK)
CODELETS_CPU(zgemm_hcore, 4, cl_zgemm_hcore_cpu_func)
// CODELETS(zgemm_hcore, 4, cl_zgemm_hcore_cpu_func, cl_zgemm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#else
CODELETS_CPU(zgemm_hcore, 7, cl_zgemm_hcore_cpu_func)
// CODELETS(zgemm_hcore, 7, cl_zgemm_hcore_cpu_func, cl_zgemm_hcore_cuda_func, STARPU_CUDA_ASYNC)
#endif
//#endif
