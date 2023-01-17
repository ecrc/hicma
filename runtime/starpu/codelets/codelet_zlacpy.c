/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dlacpy.c
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

/**
 *
 * @file codelet_zlacpy.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlacpy StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Julien Langou
 * @author Henricus Bouwmeester
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */

#include <lapacke.h>
#include <runtime/starpu/hicma_starpu.h>
#include <include/hicma_runtime_z.h>

ZCODELETS_HEADER(lacpy_hicma)

/**
 *
 * @ingroup CORE_HICMA_Complex64_t
 *
 */
void HICMA_TASK_zlacpyx(const HICMA_option_t *options,
                        HICMA_enum uplo, int m, int n, int nb,
                        int displA, const HICMA_desc_t *A, int Am, int An, int lda,
                        int displB, const HICMA_desc_t *B, int Bm, int Bn, int ldb) {
    (void) nb;
    struct starpu_codelet *codelet = &cl_zlacpy_hicma;
//    void (*callback)(void *) = options->profiling ? cl_zlacpy_callback : NULL;
    void (*callback)(void *) = NULL;

    HICMA_BEGIN_ACCESS_DECLARATION;
        HICMA_ACCESS_R(A, Am, An);
        HICMA_ACCESS_W(B, Bm, Bn);HICMA_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &uplo, sizeof(HICMA_enum),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_VALUE, &displA, sizeof(int),
            STARPU_R, RTBLKADDR(A, HICMA_Complex64_t, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_VALUE, &displB, sizeof(int),
            STARPU_W, RTBLKADDR(B, HICMA_Complex64_t, Bm, Bn),
            STARPU_VALUE, &ldb, sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zlacpy_hicma",
#endif
            0);
}

void HICMA_TASK_zlacpy(const HICMA_option_t *options,
                       HICMA_enum uplo, int m, int n, int nb,
                       const HICMA_desc_t *A, int Am, int An, int lda,
                       const HICMA_desc_t *B, int Bm, int Bn, int ldb) {
    HICMA_TASK_zlacpyx(options, uplo, m, n, nb,
                       0, A, Am, An, lda,
                       0, B, Bm, Bn, ldb);
}

#if !defined(CHAMELEON_SIMULATION)

static void cl_zlacpy_hicma_cpu_func(void *descr[], void *cl_arg) {
    HICMA_enum uplo;
    int M;
    int N;
    int displA;
    int displB;
    const HICMA_Complex64_t *A;
    int LDA;
    HICMA_Complex64_t *B;
    int LDB;

    A = (const HICMA_Complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    B = (HICMA_Complex64_t *) STARPU_MATRIX_GET_PTR(descr[1]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &displA, &LDA, &displB, &LDB);
    LAPACKE_zlacpy_work(LAPACK_COL_MAJOR, hicma_lapack_const(uplo), M, N, A + displA, LDA, B + displB, LDB);
}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlacpy_hicma, 2, cl_zlacpy_hicma_cpu_func)