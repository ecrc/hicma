/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dlaset.c
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
 * @file codelet_zlaset.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zlaset StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */

#include <lapacke.h>
#include <runtime/starpu/hicma_starpu.h>

ZCODELETS_HEADER(laset)
/**
 *
 * @ingroup CORE_HICMA_Complex64_t
 *
 *  CORE_zlaset - Sets the elements of the matrix A on the diagonal
 *  to beta and on the off-diagonals to alpha
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies which elements of the matrix are to be set
 *          = HicmaUpper: Upper part of A is set;
 *          = HicmaLower: Lower part of A is set;
 *          = HicmaUpperLower: ALL elements of A are set.
 *
 * @param[in] M
 *          The number of rows of the matrix A.  M >= 0.
 *
 * @param[in] N
 *         The number of columns of the matrix A.  N >= 0.
 *
 * @param[in] alpha
 *         The constant to which the off-diagonal elements are to be set.
 *
 * @param[in] beta
 *         The constant to which the diagonal elements are to be set.
 *
 * @param[in,out] A
 *         On entry, the M-by-N tile A.
 *         On exit, A has been set accordingly.
 *
 * @param[in] LDA
 *         The leading dimension of the array A.  LDA >= max(1,M).
 *
 */
void HICMA_TASK_zlaset(const HICMA_option_t *options,
                             HICMA_enum uplo, int M, int N,
                             HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                             const HICMA_desc_t *A, int Am, int An, int LDA) {

    struct starpu_codelet *codelet = &cl_zlaset;
    void (*callback)(void *) =  NULL;
//    void (*callback)(void *) = options->profiling ? cl_zlaset_callback : NULL;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(A, Am, An);
    HICMA_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &uplo, sizeof(HICMA_enum),
            STARPU_VALUE, &M, sizeof(int),
            STARPU_VALUE, &N, sizeof(int),
            STARPU_VALUE, &alpha, sizeof(HICMA_Complex64_t),
            STARPU_VALUE, &beta, sizeof(HICMA_Complex64_t),
            STARPU_W, RTBLKADDR(A, HICMA_Complex64_t, Am, An),
            STARPU_VALUE, &LDA, sizeof(int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zlaset_hicma",
#endif
            0);
}


#if !defined(CHAMELEON_SIMULATION)

static void cl_zlaset_cpu_func(void *descr[], void *cl_arg) {
    HICMA_enum uplo;
    int M;
    int N;
    HICMA_Complex64_t alpha;
    HICMA_Complex64_t beta;
    HICMA_Complex64_t *A;
    int LDA;

    A = (HICMA_Complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &M, &N, &alpha, &beta, &LDA);

    LAPACKE_zlaset_work(LAPACK_COL_MAJOR, hicma_lapack_const(uplo),M, N, alpha, beta, A, LDA);

}

#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zlaset,1, cl_zlaset_cpu_func)
