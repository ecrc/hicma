/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dtrsm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2019-11-21
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
 * @file codelet_dtrsm.c
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
#include <hicma_common.h>
#include <misc/auxdescutil.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(trsm_hcore_dense)

#include "flop_util_structs.h"
#include "flop_counts.h"

extern flop_counter counters[FLOP_NUMTHREADS];

#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

static int trsm_print_index_end = 0;

void
HICMA_TASK_dtrsm(const HICMA_option_t *options, HICMA_enum side, HICMA_enum uplo, HICMA_enum transA,
                       HICMA_enum diag, int m, int n, double alpha, const HICMA_desc_t *A, int Am, int An,
                       int lda, const HICMA_desc_t *B, int Bm, int Bn, int ldb) {

    struct starpu_codelet *codelet = &cl_dtrsm_hcore_dense;
    void (*callback)(void *) =  NULL;
    int sizeA = lda * m;
    int sizeB = ldb; //*nb; //@KADIR converted n to nb FIXME Size of B will be determined at runtime!!!
    int execution_rank = B->get_rankof(B, Bm, Bn);
    int rank_changed = 0;
    (void) execution_rank;

    /*  force execution on the rank owning the largest data (tile) */
    int threshold;
    char *env = getenv("HiCMA_COMM_FACTOR_THRESHOLD");
    if (env != NULL)
        threshold = (unsigned) atoi(env);
    else
        threshold = 10;
    if (sizeA > threshold * sizeB) {
        execution_rank = A->get_rankof(A, Am, An);
        rank_changed = 1;
    }
    HICMA_BEGIN_ACCESS_DECLARATION;
        HICMA_ACCESS_R(A, Am, An);
        HICMA_ACCESS_RW(B, Bm, Bn);

        if (rank_changed)
            HICMA_RANK_CHANGED(execution_rank);HICMA_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &side, sizeof(HICMA_enum),
            STARPU_VALUE, &uplo, sizeof(HICMA_enum),
            STARPU_VALUE, &transA, sizeof(HICMA_enum),
            STARPU_VALUE, &diag, sizeof(HICMA_enum),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_VALUE, &alpha, sizeof(double),
            STARPU_R, RTBLKADDR(A, double, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_RW, RTBLKADDR(B, double, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_VALUE, &Am, sizeof(int),
            STARPU_VALUE, &An, sizeof(int),
            STARPU_VALUE, &Bm, sizeof(int),
            STARPU_VALUE, &Bn, sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(HICMA_USE_MPI)
            STARPU_EXECUTE_ON_NODE, execution_rank,
#endif
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_dtrsm",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)

static void cl_dtrsm_dense_hcore_cpu_func(void *descr[], void *cl_arg) {
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
#ifdef HICMA_DISABLE_HCORE_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday(&tvalBefore, NULL);
    HICMA_enum side;
    HICMA_enum uplo;
    HICMA_enum transA;
    HICMA_enum diag;
    int m;
    int n;
    double alpha;
    double *A;
    int lda;
    double *B;
    int ldb;
    double *Brk;
    int Am;
    int An;
    int Bm;
    int Bn;

    A = (double *) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (double *) STARPU_MATRIX_GET_PTR(descr[1]);

    starpu_codelet_unpack_args(cl_arg, &side, &uplo, &transA, &diag, &m, &n, &alpha, &lda, &ldb, &Am, &An, &Bm, &Bn);

    /*CORE_dtrsm(side, uplo,*/
    /*transA, diag,*/
    /*m, n,*/
    /*alpha, A, lda,*/
    /*B, ldb);*/
    if (HICMA_get_print_index() == 1) {
        printf("%d+TRSM\t|AD(%d,%d) BV(%d,%d) m:%d lda(11):%d ldb(12):%d\n", HICMA_My_Mpi_Rank(), Am, An, Bm, Bn, m, n,
               lda, ldb);
    }

    if (HICMA_get_print_mat() == 1) {
        printf("%d\ttrsm-input A\n", __LINE__);
        _printmat(A, m, m, lda);
        printf("%d\ttrsm-input B\n", __LINE__);
        _printmat(B, m, n, ldb);
    }

    cblas_dtrsm(CblasColMajor, (CBLAS_SIDE) side, (CBLAS_UPLO) uplo, (CBLAS_TRANSPOSE) transA,
                (CBLAS_DIAG) diag, m, n, CBLAS_SADDR(alpha), A, lda, B, ldb);

    int myid = HICMA_RUNTIME_thread_rank(NULL);

    if (side == CblasLeft)
        counters[myid].trsm += flop_counts('t', m, n, 1, 0);
    else if (side == CblasRight)
        counters[myid].trsm += flop_counts('t', m, n, 2, 0);
    else
        assert(0 == "side is not CblasLeft or CblasRight");
    if (HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1 || trsm_print_index_end) {
        gettimeofday(&tvalAfter, NULL);
        printf("%d-TRSM\t|AD(%d,%d)%dx%d-%d BV(%d,%d)%dx%d-%d m:%d\t\t\t\tTRSM: %.4f\n", HICMA_My_Mpi_Rank(), Am, An, m,
               m, lda, Bm, Bn, m, n, ldb, m,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0
        );
    }
    if (HICMA_get_print_mat() == 1) {
        printf("%d\ttrsm-output\n", __LINE__);
        _printmat(B, m, n, ldb);
    }
}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dtrsm_hcore_dense, 2, cl_dtrsm_dense_hcore_cpu_func)
