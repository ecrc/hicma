/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file dgenrhs.c
 *
 *
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Rabab Alomairy
 * @date 2018-11-08
 **/
/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>

/**
 *  HICMA_dgenrhs_Tile - Generate a random matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  HICMA_dgenrhs_Tile - Generate RHS matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[out] A
 *        Store RHS values in descriptor A  
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *
 ******************************************************************************/
#include <stdio.h>
int HICMA_dgenrhs_Tile(
        HICMA_desc_t *A)
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
    HICMA_dgenrhs_Tile_Async(
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
 *  HICMA_dgenrhs_Tile_Async - Generate a RHS random matrix by tiles.
 *
 *******************************************************************************
 * @param[out] A
 *        Store RHS values in descriptor A  
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *
 ******************************************************************************/
int HICMA_dgenrhs_Tile_Async(
        HICMA_desc_t     *A,
        HICMA_sequence_t *sequence,
        HICMA_request_t  *request)
{

    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dgenrhs_Tile", "hicma not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dgenrhs_Tile", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dgenrhs_Tile", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (hicma_desc_check(A) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dgenrhs_Tile", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (hicma_min( A->m, A->n ) == 0)
        return HICMA_SUCCESS;

    hicma_pdgenrhs(
            A,
            sequence, request);

    return HICMA_SUCCESS;
}
