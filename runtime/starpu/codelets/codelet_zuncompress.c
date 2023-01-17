/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zuncompress.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/

#include <hcore_z.h>
#include <hicma.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

ZCODELETS_HEADER(uncompress_hcore)

int gemmfrk_cl_print_index = 0;
int gemmfrk_cl_print_mat = 0;
extern int print_mat;
extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern void _printmat_complex(HICMA_Complex64_t * A, int64_t m, int64_t n, int64_t ld);

/**
 *
 * CD=AU*BV'.
 * Ranks of tiles of AU are in Ark.
 * Ranks of tiles of BV are in Brk.
 * Multiplied tiles must have same rank.
 * CD is dense output.
 *
 * @ingroup CORE_double
 *
 **/
//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zuncompress(const HICMA_option_t *options,
                            HICMA_enum transA, int transB,
                            int m, int n,
                            HICMA_Complex64_t alpha,
                            const HICMA_desc_t *AUBV,
                            const HICMA_desc_t *Ark,
                            int Am, int An, int lda,
                            HICMA_Complex64_t beta,
                            const HICMA_desc_t *CD,
                            int Cm, int Cn, int ldc
                            )
{
    int nAUBV = AUBV->nb;
    struct starpu_codelet *codelet = &cl_zuncompress_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zgemmfrk_callback : NULL;*/
    void (*callback)(void*) = NULL;
    int sizeA = lda*nAUBV;
    // I converted n to k
    int sizeC = ldc*n;
    int execution_rank = CD->get_rankof( CD, Cm, Cn );
    int rank_changed=0;
    (void)execution_rank;
    
    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char* env = getenv("CHAMELEON_COMM_FACTOR_THRESHOLD");
    
    if (env != NULL)
        threshold = (unsigned)atoi(env);
    else
        threshold = 10;
    if ( sizeA > threshold*sizeC ){
        execution_rank = AUBV->get_rankof( AUBV, Am, An );
        rank_changed = 1;
    }
    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_R(AUBV, Am, An);
    HICMA_ACCESS_R(Ark, Am, An);
    HICMA_ACCESS_RW(CD, Cm, Cn);
    if (rank_changed)
        HICMA_RANK_CHANGED(execution_rank);
    HICMA_END_ACCESS_DECLARATION;
    
    /*printf("%s %d  (%d,%d) is queued to execute on rank:%d. rank_changed:%d\n", */
    /*__func__, __LINE__, Cm, Cn, execution_rank, rank_changed );*/
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,    &transA,            sizeof(HICMA_enum),
                       STARPU_VALUE,    &transB,            sizeof(HICMA_enum),
                       STARPU_VALUE,    &m,                 sizeof(int),
                       STARPU_VALUE,    &n,                 sizeof(int),
                       STARPU_VALUE,    &alpha,             sizeof(HICMA_Complex64_t),
                       STARPU_R,         RTBLKADDR(AUBV, HICMA_Complex64_t, Am, An),
                       STARPU_R,         RTBLKADDR(Ark, double, Am, An),
                       STARPU_VALUE,    &lda,               sizeof(int),
                       STARPU_VALUE,    &beta,              sizeof(HICMA_Complex64_t),
                       STARPU_RW,        RTBLKADDR(CD, HICMA_Complex64_t, Cm, Cn),
                       STARPU_VALUE,    &ldc,               sizeof(int),
                       STARPU_VALUE,    &Am,                 sizeof(int),
                       STARPU_VALUE,    &An,                 sizeof(int),
                       STARPU_VALUE,    &Cm,                 sizeof(int),
                       STARPU_VALUE,    &Cn,                 sizeof(int),
                       STARPU_VALUE,    &nAUBV,              sizeof(int),
                       STARPU_PRIORITY,  options->priority,
                       STARPU_CALLBACK,  callback,
#if defined(HICMA_USE_MPI)
                       STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zuncompress",
#endif
                       0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zuncompress_hcore_cpu_func(void *descr[], void *cl_arg)
{
    HICMA_enum transA;
    HICMA_enum transB;
    int m;
    int n;
    HICMA_Complex64_t alpha;
    HICMA_Complex64_t *AUBV;
    double *Ark;
    int lda;
    HICMA_Complex64_t beta;
    HICMA_Complex64_t *CD;
    int ldc;
    int rk;
    double acc ;
    int nAUBV;
    
    AUBV  = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double*)STARPU_MATRIX_GET_PTR(descr[1]);
    
    CD  = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    int Am, An, Cm, Cn;
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n,  &alpha, &lda, &beta, &ldc, &Am, &An, &Cm, &Cn, &nAUBV);

    int _Ark;
    if(global_always_fixed_rank == 1){
        _Ark = global_fixed_rank;
    } else {
        _Ark = Ark[0];
    }


    if(gemmfrk_cl_print_index){
        printf("+UNCOMPRESS\t|CUV(%d,%d) AUV(%d,%d)%d\n",Cm, Cn, Am, An, _Ark);
    }
    
    HICMA_Complex64_t *AU = AUBV;
   int nAU;
   if(0)printf("%s %s %d: _Ark:%d\n", __FILE__, __func__, __LINE__, _Ark);
    nAU = nAUBV/2;
    size_t nelm_AU = (size_t)lda * (size_t)nAU;
    HICMA_Complex64_t *BV = &(AUBV[nelm_AU]);
    HCORE_zuncompress(transA, transB,
                      m, n, 
                      alpha, AU, Ark, lda,
                      BV, Ark, lda,
                      beta, CD, ldc);
    
    if(gemmfrk_cl_print_index){
        printf("%d-UNCOMPRESS\t|CUV(%d,%d) AUV(%d,%d)%g mn:%d %d ldac:%d %d\n",HICMA_My_Mpi_Rank(),Cm, Cn, Am, An, Ark[0], m, n, lda, ldc);
    }
    if(gemmfrk_cl_print_mat){
        printf("InU:");zhc_printmat(AU,  m, _Ark, lda);
        printf("InV:");zhc_printmat(BV,  n, _Ark, lda);
        printf("Out:");zhc_printmat(CD,  m, n, ldc);
    }
    
}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zuncompress_hcore, 3, cl_zuncompress_hcore_cpu_func)
// CODELETS(zuncompress_hcore, 3, cl_zuncompress_hcore_cpu_func, cl_zuncompress_hcore_cuda_func, STARPU_CUDA_ASYNC)
//#endif
