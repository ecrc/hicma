/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file dhagcm.c
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Kadir Akbudak
 * @comment This file has been taken from HICMA 1.0.0
 * @date 2018-11-08
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
 * file dhagcm.c
 *
 * MORSE computational routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * date 2010-11-15
 *
 **/

#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>

/**
 *  HICMA_dhagcm_Tile - Generate a compressed matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          Which part of matrix will be generated.
 *
 * @param[in] M
 *          The number of rows of A.
 *
 * @param[in] N
 *          The order of the matrix A. N >= 0.
 *
 * @param[out] AUV
 *          On exit, U and V factors.
 *
 * @param[out] AD
 *          On exit, diagonal tiles.
 *
 * @param[out] Ark
 *          On exit, rank for each tile.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] seed
 *          The seed used in the random generation.
 *
 * @param[in] maxrank
 *          The limit for ranks. 
 *
 * @param[in] tol
 *          Threshold used in approximation.
 *

 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_dhagcm_Tile(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank,
        double tol
) {
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dhagcm_Tile", "hicma not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);
    HICMA_dhagcm_Tile_Async(
            uplo,
            AUV,
            Ark,
            numrows_matrix,
            numcols_matrix,
            numrows_block,
            numcols_block,
            maxrank, tol,
            sequence, &request);
    hicma_sequence_wait(hicma, sequence);
    //RUNTIME_desc_getoncpu(AD);
    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/**
 *
 *  HICMA_dhagcm_Tile_Async - Generate a compressed matrix by tiles.
 *  Non-blocking equivalent of HICMA_dhagcm_Tile().
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
 ******************************************************************************/
int HICMA_dhagcm_Tile_Async(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank, double tol,
        HICMA_sequence_t *sequence,
        HICMA_request_t *request) {
    HICMA_desc_t *A = AUV; // FIXME

    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dhagcm_Tile", "hicma not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dhagcm_Tile", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dhagcm_Tile", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (hicma_desc_check(A) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dhagcm_Tile", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb) {
        hicma_error("HiCMA_dhagcm_Tile", "only square tiles supported");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }*/

    /* Quick return */
    if (hicma_min(A->m, A->n) == 0)
        return HICMA_SUCCESS;

    hicma_pdhagcm(
            uplo,
            AUV,
            Ark,
            numrows_matrix,
            numcols_matrix,
            numrows_block,
            numcols_block,
            maxrank, tol,
            sequence, request);

    return HICMA_SUCCESS;
}
