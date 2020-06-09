/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file zpotrf.c
 *
 * This file contains top-level functions for Cholesky factorization.
 *  
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 */
/*
 * file zpotrf.c
 *
 * MORSE computational routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * author Jakub Kurzak
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * date 2010-11-15
 *
 */
#include "morse.h"
#include "control/common.h"
#include "control/hicma_common.h"

int HICMA_zgene_Tile( int ntrian, int nip, int local_nt, int NB, MORSE_desc_t *A)
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_zpotrf_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
  //printf("\nHERERE1ntrian:%d, nip:%d, local_nt:%d, NB:%d\n", ntrian, nip, local_nt, NB);
    HICMA_zgene_Tile_Async( ntrian,  nip, local_nt,  NB, A, sequence, &request);

    MORSE_Desc_Flush( A, sequence );

    morse_sequence_wait(morse, sequence);

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}
/***************************************************************************//**
 *
 *  HICMA_zpotrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  positive definite positive definite matrix.
 *  Non-blocking equivalent of HICMA_zpotrf_Tile().
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
int HICMA_zgene_Tile_Async( int ntrian, int nip, int local_nt, int NB, MORSE_desc_t *A, MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_zgene_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("HICMA_zgene_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("HICMA_zgene_Tile_Async", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (
            (morse_desc_check(A) != MORSE_SUCCESS)
            ){
        morse_error("HICMA_zgene_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Quick return */
/*
    if (chameleon_max(N, 0) == 0)
        return MORSE_SUCCESS;
*/
   // printf("\nHERERE2ntrian:%d, nip:%d, local_nt:%d, NB:%d\n", ntrian, nip, local_nt, NB);
    compute_gene( ntrian, nip, local_nt, NB, A,
                    sequence, &request );

    return MORSE_SUCCESS;
}
