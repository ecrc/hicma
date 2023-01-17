/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file zplrnt.c
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
 * file zplrnt.c
 *
 * @brief Chameleon zplrnt parallel algorithm
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
#include <control/hicma_compute_z.h>
#include <include/hicma_z.h>

/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t
 *
 *  HICMA_zplrnt - Generate a random matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *******************************************************************************
 *
 * @sa HICMA_zplrnt_Tile
 * @sa HICMA_zplrnt_Tile_Async
 * @sa HICMA_cplrnt
 * @sa HICMA_dplrnt
 * @sa HICMA_splrnt
 * @sa HICMA_zplghe
 * @sa HICMA_zplgsy
 *
 */
int HICMA_zplrnt(int M, int N,
                 HICMA_Complex64_t *A, int LDA,
                 unsigned long long int seed) {
    int NB;
    int status;
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    HICMA_desc_t descAl, descAt;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zplrnt", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    /* Check input arguments */
    if (M < 0) {
        hicma_error("HiCMA_zplrnt", "illegal value of M");
        return -1;
    }
    if (N < 0) {
        hicma_error("HiCMA_zplrnt", "illegal value of N");
        return -2;
    }
    if (LDA < hicma_max(1, M)) {
        hicma_error("HiCMA_zplrnt", "illegal value of LDA");
        return -4;
    }
    /* Quick return */
    if (hicma_min(M, N) == 0)
        return HICMA_SUCCESS;

    /* Tune NB depending on M, N & NRHS; Set NBNB */
    status = hicma_tune(HICMA_FUNC_ZGEMM, M, N, 0);
    if (status != HICMA_SUCCESS) {
        hicma_error("HiCMA_zplrnt", "hicma_tune() failed");
        return status;
    }

    /* Set NT */
    NB = HICMA_NB;
    hicma_sequence_create(hicma, &sequence);

    /* Submit the matrix conversion */
    hicma_zlap2tile(hicma, &descAl, &descAt, HicmaDescOutput, HicmaUpperLower,
                    A, NB, NB, LDA, N, M, N, sequence, &request);

    /* Call the tile interface */
    HICMA_zplrnt_Tile_Async(&descAt, seed, sequence, &request);

    /* Submit the matrix conversion back */
    hicma_ztile2lap(hicma, &descAl, &descAt,
                    HicmaDescOutput, HicmaUpperLower, sequence, &request);

    hicma_sequence_wait(hicma, sequence);

    /* Cleanup the temporary data */
    hicma_ztile2lap_cleanup(hicma, &descAl, &descAt);

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t_Tile
 *
 *  HICMA_zplrnt_Tile - Generate a random matrix by tiles.
 *  Tile equivalent of HICMA_zplrnt().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          On exit, The random matrix A generated.
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa HICMA_zplrnt
 * @sa HICMA_zplrnt_Tile_Async
 * @sa HICMA_cplrnt_Tile
 * @sa HICMA_dplrnt_Tile
 * @sa HICMA_splrnt_Tile
 * @sa HICMA_zplghe_Tile
 * @sa HICMA_zplgsy_Tile
 *
 */
int HICMA_zplrnt_Tile(HICMA_desc_t *A,
                      unsigned long long int seed) {
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zplrnt_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);

    HICMA_zplrnt_Tile_Async(A, seed, sequence, &request);

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
 *  HICMA_zplrnt_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of HICMA_zplrnt_Tile().
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
 * @sa HICMA_zplrnt
 * @sa HICMA_zplrnt_Tile
 * @sa HICMA_cplrnt_Tile_Async
 * @sa HICMA_dplrnt_Tile_Async
 * @sa HICMA_splrnt_Tile_Async
 * @sa HICMA_zplghe_Tile_Async
 * @sa HICMA_zplgsy_Tile_Async
 *
 */
int HICMA_zplrnt_Tile_Async(HICMA_desc_t *A,
                            unsigned long long int seed,
                            HICMA_sequence_t *sequence,
                            HICMA_request_t *request) {
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zplrnt_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_zplrnt_Tile", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_zplrnt_Tile", "NULL request");
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
        hicma_error("HiCMA_zplrnt_Tile", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (A->nb != A->mb) {
        hicma_error("HiCMA_zplrnt_Tile", "only square tiles supported");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (hicma_min(A->m, A->n) == 0)
        return HICMA_SUCCESS;

    hicma_pzplrnt(A, seed, sequence, request);

    return HICMA_SUCCESS;
}
