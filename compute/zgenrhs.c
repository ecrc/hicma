/*
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file zgenrhs.c
 *
 *
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Rabab Alomairy
 * @date 2018-11-08
 **/
/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */
#include "morse.h"
#include "control/common.h"
#include "control/hicma_common.h"

/**
 *  HICMA_zgenrhs_Tile - Generate a random matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  HICMA_zgenrhs_Tile - Generate RHS matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[out] A
 *        Store RHS values in descriptor A  
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *
 ******************************************************************************/
#include <stdio.h>
int HICMA_zgenrhs_Tile(
        MORSE_desc_t *A)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgytlr_Tile", "morse not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    HICMA_zgenrhs_Tile_Async(
            A,
            sequence, &request );
    MORSE_Desc_Flush( A, sequence );
    morse_sequence_wait(morse, sequence);
    //RUNTIME_desc_getoncpu(AD);

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);

    return status;
}

/**
 *
 *  HICMA_zgenrhs_Tile_Async - Generate a RHS random matrix by tiles.
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
int HICMA_zgenrhs_Tile_Async(
        MORSE_desc_t     *A,
        MORSE_sequence_t *sequence,
        MORSE_request_t  *request)
{

    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgenrhs_Tile", "morse not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgenrhs_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgenrhs_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zgenrhs_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return MORSE_SUCCESS;

    hicma_pzgenrhs(
            A,
            sequence, request);

    return MORSE_SUCCESS;
}
