/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dgemm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Ali Charara
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

/*
 *
 * @file codelet_dgemm.c
 *
 *  MORSE computational routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-08
 * @precisions normal z -> s d c
 *
 **/

#include <sys/time.h>
#include <hicma.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(gemm_hcore_dense)

#include "hcore_d.h"

extern flop_counter counters[FLOP_NUMTHREADS];

extern int global_always_fixed_rank;
extern int global_fixed_rank;
extern int print_index;
extern int print_mat;

extern void _printmat(double *A, int64_t m, int64_t n, int64_t ld);

/**
 *
 * @ingroup hcore_dgemm
 *
 **/

void HICMA_TASK_dgemm(const HICMA_option_t *options,
                            HICMA_enum transA, int transB,
                            int m, int n, int k,
                            double alpha,
                            const HICMA_desc_t *A,
                            int Am, int An, int lda,
                            const HICMA_desc_t *B,
                            int Bm, int Bn, int ldb,
                            double beta,
                            const HICMA_desc_t *C,
                            int Cm, int Cn, int ldc) {
    int nA = A->nb;
    int nB = B->nb;
    int nC = C->nb;
    struct starpu_codelet *codelet = &cl_dgemm_hcore_dense;
    /*void (*callback)(void*) = options->profiling ? cl_dgemm_hcore_callback : NULL;*/
    void (*callback)(void *) =  NULL;

    int sizeA = lda * nA; //FIXME Think about scheduling of tasks according to sizes of the matrices
    int sizeB = ldb * nB;
    int sizeC = ldc * nC;
    int execution_rank = C->get_rankof(C, Cm, Cn);
    int rank_changed = 0;
    (void) execution_rank;

    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char *env = getenv("HiCMA_COMM_FACTOR_THRESHOLD");

    int ifval = 0, elseifval = 0, initialval = execution_rank;
    if (env != NULL)
        threshold = (unsigned) atoi(env);
    else
        threshold = 10;
    if (sizeA > threshold * sizeC) {
        execution_rank = A->get_rankof(A, Am, An);
        ifval = execution_rank;
        rank_changed = 1;
    } else if (sizeB > threshold * sizeC) {
        execution_rank = B->get_rankof(B, Bm, Bn);
        elseifval = execution_rank;
        rank_changed = 1;
    }

    //printf("%d,%d %d,%d %d,%d\n", Am, An, Bm, Bn, Cm, Cn);
    //printf("initialval:\t%d if:%d\t else:\t%d rc:\t%d\n", initialval, ifval, elseifval, rank_changed);
    HICMA_BEGIN_ACCESS_DECLARATION;
        HICMA_ACCESS_R(A, Am, An);
        HICMA_ACCESS_R(B, Bm, Bn);
        HICMA_ACCESS_RW(C, Cm, Cn);

        if (rank_changed)
            HICMA_RANK_CHANGED(execution_rank);HICMA_END_ACCESS_DECLARATION;

    //printf("%s %d n:%d\n", __func__, __LINE__,n );
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &transA, sizeof(HICMA_enum),
            STARPU_VALUE, &transB, sizeof(HICMA_enum),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_VALUE, &k, sizeof(int),
            STARPU_VALUE, &alpha, sizeof(double),
            STARPU_R, RTBLKADDR(A, double, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_R, RTBLKADDR(B, double, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_VALUE, &beta, sizeof(double),
            STARPU_RW, RTBLKADDR(C, double, Cm, Cn),
            STARPU_VALUE, &ldc, sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(HICMA_USE_MPI)
            STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_dgemm",
#endif
            0);
}

#if !defined(CHAMELEON_SIMULATION)

static void cl_dgemm_dense_hcore_cpu_func(void *descr[], void *cl_arg) {
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday(&tvalBefore, NULL);
    HICMA_enum transA;
    HICMA_enum transB;
    int m;
    int n;
    int k;
    double alpha;
    double *A = NULL;
    int lda;
    double *B = NULL;
    int ldb;
    double beta;
    double *C = NULL;
    int ldc;

    int idescr = 0;
    A = (double *) STARPU_MATRIX_GET_PTR(descr[idescr++]);
    B = (double *) STARPU_MATRIX_GET_PTR(descr[idescr++]);
    C = (double *) STARPU_MATRIX_GET_PTR(descr[idescr++]);

    starpu_codelet_unpack_args(cl_arg, &transA, &transB, &m, &n, &k, &alpha, &lda, &ldb, &beta, &ldc);

    //printf("%d,%d:%g\t%d,%d:%g\t%d,%d:%g\n", Am, An, *Ark, Bm, Bn, *Brk, Cm, Cn, *Crk);

    char datebuf_start[128];
    datebuf_start[0] = '\0';
    if (print_index) {
        time_t timer;
        struct tm *tm_info;
        gettimeofday(&tvalAfter, NULL);
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf_start, 26, "%Y-%m-%d %H:%M:%S", tm_info); \
//        printf("%d+GEMM\t|CUV(%d,%d) AUV(%d,%d) BUV(%d,%d)\t\t\t\t\tGEMM: %s\n", HICMA_My_Mpi_Rank(), Cm, Cn,
//               Am, An, Bm, Bn, datebuf_start);
    }

    flop_counter flops;
    flops.update = 0;

    HCORE_dgemm_dense((HICMA_enum) transA, (HICMA_enum) transB, m, n, k, alpha, A, lda, B, ldb, beta, C, ldc);

    int myid = HICMA_RUNTIME_thread_rank(NULL);
    counters[myid].update += flops.update;

    if (print_index) {
        char datebuf[128];
        time_t timer;
        struct tm *tm_info;
        gettimeofday(&tvalAfter, NULL);
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S", tm_info); \
//        printf("%d-GEMM\t|CUV(%d,%d) AUV(%d,%d) BUV(%d,%d) acc:%e rk:%d maxrk:%d\t\t\tGEMM: %.4f\t%s---%s\n",
//               HICMA_My_Mpi_Rank(), Cm, Cn, Am, An, Bm, Bn,
//               (tvalAfter.tv_sec - tvalBefore.tv_sec)
//               + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0,
//               datebuf_start, datebuf
//        );
    }
}

#endif /* !defined(HICMA_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dgemm_hcore_dense, 3, cl_dgemm_dense_hcore_cpu_func)
