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

/***************************************************************************//**
 *
 *  HICMA_zpotrf_Tile - Computes the Cholesky factorization of a symmetric 
 *  positive definite matrix in tile low-rank (TLR) format.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = MorseUpper: Upper triangle of A is stored (Not supported yet)
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite TLR matrix A.
 *          If uplo = MorseUpper, the leading N-by-N upper triangular part of A
 *          contains the upper triangular part of the matrix A, and the strictly lower triangular
 *          part of A is not referenced.
 *          If UPLO = 'L', the leading N-by-N lower triangular part of A contains the lower
 *          triangular part of the matrix A, and the strictly upper triangular part of A is not
 *          referenced.
 *          On exit, if return value = 0, the factor U or L from the Cholesky factorization
 *          A = L*L**H.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 ******************************************************************************/
int HICMA_zpotrf_Tile(MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        int rk, int maxrk, double acc
        )
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
    HICMA_zpotrf_Tile_Async(uplo,
            AUV, AD, Ark,
            rk, maxrk, acc,
            sequence, &request
            );
    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

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
int HICMA_zpotrf_Tile_Async(MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        int rk, int maxrk, double acc,
        MORSE_sequence_t *sequence, MORSE_request_t *request
        )
{
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_zpotrf_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("HICMA_zpotrf_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("HICMA_zpotrf_Tile_Async", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (
            (morse_desc_check(AUV) != MORSE_SUCCESS)
            || (morse_desc_check(AD) != MORSE_SUCCESS)
            || (morse_desc_check(Ark) != MORSE_SUCCESS)
            ){
        morse_error("HICMA_zpotrf_Tile_Async", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (AD->nb != AD->mb) {
        morse_error("HICMA_zpotrf_Tile_Async", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (uplo != MorseUpper && uplo != MorseLower) {
        morse_error("HICMA_zpotrf_Tile_Async", "illegal value of uplo");
        return morse_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (chameleon_max(N, 0) == 0)
        return MORSE_SUCCESS;
*/

    hicma_pzpotrf(uplo, AUV, AD, Ark, sequence, request,
            rk, maxrk, acc
            );

    return MORSE_SUCCESS;
}
