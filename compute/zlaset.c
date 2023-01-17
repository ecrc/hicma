/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file zlaset.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 **/

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

/**
 *
 * file zlaset.c
 *
 * @brief Chameleon zlaset parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/

#include <control/common.h>
#include <include/hicma_z.h>
#include <control/hicma_compute_z.h>
/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t
 *
 *  HICMA_zlaset copies all or part of a two-dimensional matrix A to another
 *  matrix B
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = HicmaUpperLower: All the matrix A
 *            = HicmaUpper: Upper triangular part is set. The lower
 *            triangle is unchanged.
 *            = HicmaLower: Lower triangular part is set. The upper
 *            triangle is unchange.
 *
 * @param[in] M
 *          The number of rows of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of columns of the matrix A. N >= 0.
 *
 * @param[in] alpha
 *          All the offdiagonal array elements are set to alpha.
 *
 * @param[in] beta
 *          All the diagonal array elements are set to beta.
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 *******************************************************************************
 *
 * @sa HICMA_zlaset_Tile
 * @sa HICMA_zlaset_Tile_Async
 * @sa HICMA_claset
 * @sa HICMA_dlaset
 * @sa HICMA_slaset
 *
 */
int HICMA_zlaset(HICMA_enum uplo, int M, int N,
                       HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                       HICMA_Complex64_t *A, int LDA) {
    int NB;
    int status;
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    HICMA_desc_t descAl, descAt;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zlaset", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if ((uplo != HicmaUpperLower) &&
        (uplo != HicmaUpper) &&
        (uplo != HicmaLower)) {
        hicma_error("HiCMA_zlaset", "illegal value of uplo");
        return -1;
    }
    if (M < 0) {
        hicma_error("HiCMA_zlaset", "illegal value of M");
        return -2;
    }
    if (N < 0) {
        hicma_error("HiCMA_zlaset", "illegal value of N");
        return -3;
    }
    if (LDA < hicma_max(1, M)) {
        hicma_error("HiCMA_zlaset", "illegal value of LDA");
        return -5;
    }

    /* Quick return */
    if (hicma_min(N, M) == 0)
        return (double) 0.0;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = hicma_tune(HICMA_FUNC_ZGEMM, M, N, 0);
    if (status != HICMA_SUCCESS) {
        hicma_error("HiCMA_zlaset", "hicma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = HICMA_NB;

    hicma_sequence_create(hicma, &sequence);

    /* Submit the matrix conversion */
    hicma_zlap2tile(hicma, &descAl, &descAt, HicmaDescInout, uplo,
                    A, NB, NB, LDA, N, M, N, sequence, &request);

    /* Call the tile interface */
    HICMA_zlaset_Tile_Async(uplo, alpha, beta, &descAt, sequence, &request);

    /* Submit the matrix conversion back */
    hicma_ztile2lap(hicma, &descAl, &descAt,
                    HicmaDescInout, uplo, sequence, &request);

    hicma_sequence_wait(hicma, sequence);

    /* Cleanup the temporary data */
    hicma_ztile2lap_cleanup(hicma, &descAl, &descAt);

    hicma_sequence_destroy(hicma, sequence);
    return HICMA_SUCCESS;
}

/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t_Tile
 *
 *  HICMA_zlaset_Tile - Tile equivalent of HICMA_zlaset().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Specifies the part of the matrix A to be copied to B.
 *            = HicmaUpperLower: All the matrix A
 *            = HicmaUpper: Upper triangular part
 *            = HicmaLower: Lower triangular part
 *
 * @param[in,out] A
 *          On entry, the m by n matrix A.
 *          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
 *                   A(i,i) = BETA,  1 <= i <= min(m,n)
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa HICMA_zlaset
 * @sa HICMA_zlaset_Tile_Async
 * @sa HICMA_claset_Tile
 * @sa HICMA_dlaset_Tile
 * @sa HICMA_slaset_Tile
 *
 */
int HICMA_zlaset_Tile(HICMA_enum uplo,
                            HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                            HICMA_desc_t *A) {
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zlaset_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);

    HICMA_zlaset_Tile_Async(uplo, alpha, beta, A, sequence, &request);

    HICMA_Desc_Flush(A, sequence);

    hicma_sequence_wait(hicma, sequence);
    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t_Tile_Async
 *
 *  HICMA_zlaset_Tile_Async - Non-blocking equivalent of HICMA_zlaset_Tile().
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 *******************************************************************************
 *
 * @sa HICMA_zlaset
 * @sa HICMA_zlaset_Tile
 * @sa HICMA_claset_Tile_Async
 * @sa HICMA_dlaset_Tile_Async
 * @sa HICMA_slaset_Tile_Async
 *
 */
int HICMA_zlaset_Tile_Async(HICMA_enum uplo,
                                  HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                                  HICMA_desc_t *A,
                                  HICMA_sequence_t *sequence, HICMA_request_t *request) {
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zlaset_Tile_Async", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_zlaset_Tile_Async", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_zlaset_Tile_Async", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS) {
        request->status = HICMA_SUCCESS;
    } else {
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);
    }

    /* Check descriptors for correctness */
    if (hicma_desc_check(A) != HICMA_SUCCESS) {
        hicma_error("HiCMA_zlaset_Tile_Async", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        hicma_error("HiCMA_zlaset_Tile_Async", "only square tiles supported");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((uplo != HicmaUpperLower) &&
        (uplo != HicmaUpper) &&
        (uplo != HicmaLower)) {
        hicma_error("HiCMA_zlaset_Tile_Async", "illegal value of uplo");
        return -1;
    }
    /* Quick return */
    if (hicma_min(A->m, A->n) == 0) {
        return HICMA_SUCCESS;
    }

    hicma_pzlaset(uplo, alpha, beta, A, sequence, request);

    return HICMA_SUCCESS;
}
