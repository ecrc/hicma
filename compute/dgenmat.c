/*
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file dgenmat.c
 *
 * This file contains tile low-rank (TLR) matrix generation functions.
 *
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Rabab Alomairy
 * @author Kadir Akbudak
 * @date 2020-03-04
 **/

#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>

/**
 *  HICMA_dgenmat_Tile - Generate application matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *
 *******************************************************************************
 *
 * @param[out] A
 *        Store application values in descriptor A  
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *
 ******************************************************************************/
#include <stdio.h>
int HICMA_dgenmat_Tile(
        HICMA_desc_t *A
        )
{
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dgytlr_Tile", "hicma not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);
    HICMA_dgenmat_Tile_Async(
            A,
            sequence, &request );
    HICMA_Desc_Flush( A, sequence );
    hicma_sequence_wait(hicma, sequence);
    //RUNTIME_desc_getoncpu(AD);

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);

    return status;
}

/**
 *
 *  HICMA_dgenmat_Tile_Async - Generate application matrix by tiles.
 *  May return before the computation is finished.
 *  Allows for pipelining of operations at runtime.
 *
 *******************************************************************************
 *
 * @param[out] A
 *        Store application values in descriptor A  
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/
int HICMA_dgenmat_Tile_Async(
        HICMA_desc_t *A,
        HICMA_sequence_t *sequence,
        HICMA_request_t  *request)
{

    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dgenmat_Tile_Async", "hicma not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dgenmat_Tile_Async", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dgenmat_Tile_Async", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (hicma_desc_check(A) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dgytlr_Tile", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb) {
        hicma_error("HiCMA_dgytlr_Tile", "only square tiles supported");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }*/

    /* Quick return */
    if (hicma_min( A->m, A->n ) == 0)
        return HICMA_SUCCESS;

    hicma_pdgenmat(
            A,
            sequence, request);

    return HICMA_SUCCESS;
}
