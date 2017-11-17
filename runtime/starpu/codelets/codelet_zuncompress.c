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
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
/*#include "runtime/starpu/include/runtime_codelet_z.h"*/

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(uncompress_hcore)

#include "hcore_z.h"

int gemmfrk_cl_print_index = 0;
int gemmfrk_cl_print_mat = 0;
extern int print_mat;
extern void _printmat(double * A, int64_t m, int64_t n, int64_t ld);
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

void HICMA_TASK_zuncompress(const MORSE_option_t *options,
        MORSE_enum transA, int transB,
        int m, int n, 
        double alpha,
        const MORSE_desc_t *AUBV,
        const MORSE_desc_t *Ark,
        int Am, int An, int lda,
        double beta,
        const MORSE_desc_t *CD,
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
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(AUBV, Am, An);
    MORSE_ACCESS_R(Ark, Am, An);
    MORSE_ACCESS_RW(CD, Cm, Cn);
    if (rank_changed)
        MORSE_RANK_CHANGED(execution_rank);
    MORSE_END_ACCESS_DECLARATION;

    /*printf("%s %d  (%d,%d) is queued to execute on rank:%d. rank_changed:%d\n", */
    /*__func__, __LINE__, Cm, Cn, execution_rank, rank_changed );*/
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &transA,            sizeof(MORSE_enum),
            STARPU_VALUE,    &transB,            sizeof(MORSE_enum),
            STARPU_VALUE,    &m,                 sizeof(int),
            STARPU_VALUE,    &n,                 sizeof(int),
            STARPU_VALUE,    &alpha,             sizeof(double),
            STARPU_R,         RTBLKADDR(AUBV, double, Am, An),
            STARPU_R,         RTBLKADDR(Ark, double, Am, An),
            STARPU_VALUE,    &lda,               sizeof(int),
            STARPU_VALUE,    &beta,              sizeof(double),
            STARPU_RW,        RTBLKADDR(CD, double, Cm, Cn),
            STARPU_VALUE,    &ldc,               sizeof(int),
            STARPU_VALUE,    &Am,                 sizeof(int),
            STARPU_VALUE,    &An,                 sizeof(int),
            STARPU_VALUE,    &Cm,                 sizeof(int),
            STARPU_VALUE,    &Cn,                 sizeof(int),
            STARPU_VALUE,    &nAUBV,              sizeof(int),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
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
    MORSE_enum transA;
    MORSE_enum transB;
    int m;
    int n;
    double alpha;
    double *AUBV;
    double *Ark;
    int lda;
    double beta;
    double *CD;
    int ldc;
    int rk;
    double acc ;
    int nAUBV;

    AUBV  = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[1]);

    CD  = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    int Am, An, Bm, Bn, Cm, Cn;
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n,  &alpha, &lda, &beta, &ldc, &Am, &An, &Cm, &Cn, &nAUBV);
    if(gemmfrk_cl_print_index){
        //printf("+GEMMFRK\t|CUV(%d,%d) AUV(%d,%d) BUV(%d,%d)\n",Cm, Cn, Am, An, Bm, Bn);
    }
    if(gemmfrk_cl_print_index){
        printf("%d-UNCOMPRESS\t|CUV(%d,%d) AUV(%d,%d)%g mn:%d %d ldac:%d %d\n",MORSE_My_Mpi_Rank(),Cm, Cn, Am, An, Ark[0], m, n, lda, ldc);
    }

    double *AU = AUBV;
    int nAU = nAUBV/2;
    size_t nelm_AU = (size_t)lda * (size_t)nAU;
    double *BV = &(AUBV[nelm_AU]);
    HCORE_zuncompress(transA, transB,
            m, n, 
            alpha, AU, Ark, lda,
             BV, Ark, lda,
            beta, CD, ldc);


}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zuncompress_hcore, 3, cl_zuncompress_hcore_cpu_func)
// CODELETS(zuncompress_hcore, 3, cl_zuncompress_hcore_cpu_func, cl_zuncompress_hcore_cuda_func, STARPU_CUDA_ASYNC)
