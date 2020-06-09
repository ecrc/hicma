/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zgemm_bdcd.c
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
ZCODELETS_HEADER(gemmbdcd_hcore)

//UPDATE this definition. I only copy-paste from runtime/starpu/codelets/codelet_zcallback.c
/*CHAMELEON_CL_CB(zgemm_hcore,         starpu_matrix_get_nx(task->handles[2]), starpu_matrix_get_ny(task->handles[2]), starpu_matrix_get_ny(task->handles[0]),     2. *M*N*K) [> If A^t, computation is wrong <]*/

#include "hcore_z.h"

extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern int print_mat;
extern void _printmat(double * A, int64_t m, int64_t n, int64_t ld);
/**
 *
 * @ingroup hcore_zgemmbdcd
 *
 **/

void HICMA_TASK_zgemm_bdcd(const MORSE_option_t *options,
        MORSE_enum transA, int transB,
        int m, int n,
        double alpha,
        const MORSE_desc_t *AUV,
        const MORSE_desc_t *Ark,
        int Am, int An, int lda,
        const MORSE_desc_t *BD,
        int Bm, int Bn, int ldb,
        double beta,
        const MORSE_desc_t *CD,
        int Cm, int Cn, int ldc
        )
{
    int nAUV = AUV->nb;
    struct starpu_codelet *codelet = &cl_zgemmbdcd_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zgemm_hcore_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    MORSE_starpu_ws_t *h_work = (MORSE_starpu_ws_t*)(options->ws_host);
    /*printf("%s %d:\t%p %p\n", __FILE__, __LINE__, h_work, options->ws_host);*/

    int sizeA = lda*nAUV; //FIXME Think about scheduling of tasks according to sizes of the matrices
    int sizeB = ldb*n;
    int sizeC = ldc*m; 
    int execution_rank = CD->get_rankof( CD, Cm, Cn );
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
        execution_rank = BD->get_rankof( BD, Bm, Bn );
        elseifval = execution_rank;
        rank_changed = 1;
    }
    //printf("m:%d n:%d k:%d nb:%d\n", m, n, k, nb); all of them are nb (1156)
    //printf("initialval:\t%d if:%d\t else:\t%d rc:\t%d\n", initialval, ifval, elseifval, rank_changed);
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_R(AUV, Am, An);
    MORSE_ACCESS_R(BD, Bm, Bn);
    MORSE_ACCESS_RW(CD, Cm, Cn);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    MORSE_ACCESS_R(Ark, Am, An);
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
            STARPU_VALUE,     &alpha,             sizeof(double),
            STARPU_R,         RTBLKADDR(AUV, double, Am, An),
            STARPU_VALUE,     &lda,               sizeof(int),
            STARPU_R,         RTBLKADDR(BD, double, Bm, Bn),
            STARPU_VALUE,     &ldb,               sizeof(int),
            STARPU_VALUE,     &beta,              sizeof(double),
            STARPU_RW,        RTBLKADDR(CD, double, Cm, Cn),
#if !defined(HICMA_ALWAYS_FIX_RANK)
            STARPU_R,         RTBLKADDR(Ark, double, Am, An),
#endif
            STARPU_VALUE,     &ldc,               sizeof(int),
            STARPU_VALUE,     &Am,                 sizeof(int),
            STARPU_VALUE,     &An,                 sizeof(int),
            STARPU_VALUE,     &Bm,                 sizeof(int),
            STARPU_VALUE,     &Bn,                 sizeof(int),
            STARPU_VALUE,     &Cm,                 sizeof(int),
            STARPU_VALUE,     &Cn,                 sizeof(int),
            STARPU_VALUE,     &nAUV,                 sizeof(int),
            STARPU_SCRATCH,   options->ws_worker,
            STARPU_VALUE,     &h_work,            sizeof(MORSE_starpu_ws_t *),
            STARPU_PRIORITY,  options->priority,
            STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_USE_MPI)
            STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_zgemm_bdcd",
#endif
            0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgemmbdcd_hcore_cpu_func(void *descr[], void *cl_arg)
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
    double alpha;
    double *AUV = NULL;
    double *Ark = NULL;
    int lda;
    double *BD = NULL;
    int ldb;
    double beta;
    double *CD = NULL;
    int ldc;
    int nAUV;

    int idescr = 0;
    AUV  = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    BD  = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
    CD  = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#if !defined(HICMA_ALWAYS_FIX_RANK)
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);
#else
    double _gemm_rank = global_fixed_rank;
    Ark = &_gemm_rank;
#endif

    double* work = NULL;
    work = (double *)STARPU_MATRIX_GET_PTR(descr[idescr++]);

    int Am, An, Bm, Bn, Cm, Cn;

    MORSE_starpu_ws_t *h_work;
    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &alpha, &lda, &ldb, &beta, &ldc, &Am, &An, &Bm, &Bn, &Cm, &Cn, &nAUV, &h_work);

    double *AU = AUV;

    int nAU = nAUV/2;
    size_t nelm_AU = (size_t)lda * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);

    char datebuf_start[128];
    if(HICMA_get_print_index()){
        time_t timer;
        struct tm* tm_info;
        gettimeofday (&tvalAfter, NULL);
        time(&timer); 
        tm_info = localtime(&timer); 
        strftime(datebuf_start, 26, "%Y-%m-%d %H:%M:%S",tm_info); 
        printf("%d+GEMMBDCD\t|CD(%d,%d) AUV(%d,%d)%g BD(%d,%d) m:%d n:%d lda:%d ldb:%d ldc:%d \t\t\t\t\tGEMMBDCD: %s\n",MORSE_My_Mpi_Rank(), Cm, Cn, Am, An, Ark[0], Bm, Bn, m, n, lda, ldb, ldc, datebuf_start);
    }

    int isTransA = transA == MorseTrans;
    int isTransB = transB == MorseTrans;
    if(isTransB == 1){
        printf("%s %d %s: Transpose of B is not supported yet. isTransB: %d transB:%d\n", __FILE__, __LINE__, __func__, isTransB, transB);
        exit(101);
    }

    HCORE_zgemmbdcd(transA, transB,
            m, n,
            alpha, (isTransA ? AV : AU), (isTransA ? AU : AV), Ark, lda,
            BD, ldb,
            beta, CD, ldc, work);

    if(HICMA_get_print_index() || HICMA_get_print_index_end()){
        char datebuf[128];
        time_t timer;
        struct tm* tm_info;
        gettimeofday (&tvalAfter, NULL);
        time(&timer); 
        tm_info = localtime(&timer); 
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S",tm_info); 
        printf("%d-GEMMBDCD\t|CD(%d,%d) AUV(%d,%d)%g BD(%d,%d)\t\t\tGEMMBDCD: %.4f\t%s---%s\n",MORSE_My_Mpi_Rank(),Cm, Cn, Am, An, Ark[0], Bm, Bn,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0,
                 datebuf_start, datebuf
                );
    }
}
#endif /* !defined(MORSE_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgemmbdcd_hcore, 5, cl_zgemmbdcd_hcore_cpu_func)

