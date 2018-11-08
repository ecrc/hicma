/*
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file zgytlr.c
 *
 * This file contains tile low-rank (TLR) matrix generation functions.
 *
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
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
 */
/*
 * file zplrnt.c
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
 * precisions normal z -> s d c
 *
 */
#include "morse.h"
#include "control/common.h"
#include "control/hicma_common.h"

/**
 *  HICMA_zgytlr_Tile - Generate a random matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  HICMA_zgytlr_Tile - Generate matrix by tiles.
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
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *          \retval <0 if -i, the i-th argument had an illegal value
 *
 *
 ******************************************************************************/
#include <stdio.h>
int HICMA_zgytlr_Tile(
        MORSE_enum   uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        unsigned long long int seed,
        int maxrank,
        double tol,
        int compress_diag,
        MORSE_desc_t *Dense
        )
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
    HICMA_zgytlr_Tile_Async(
            uplo,
            AUV,
            AD,
            Ark,
            seed, maxrank, tol,
            compress_diag,
            Dense,
            sequence, &request );
    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( AUV, sequence ); //added due to stall when parsec is used
    MORSE_Desc_Flush( Ark, sequence ); //added due to stall when parsec is used
    MORSE_Desc_Flush( Dense, sequence ); //added due to stall when parsec is used
    morse_sequence_wait(morse, sequence);
    //RUNTIME_desc_getoncpu(AD);

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);

    return status;
}

/**
 *
 *  HICMA_zgytlr_Tile_Async - Generate a random matrix by tiles.
 *  Non-blocking equivalent of HICMA_zgytlr_Tile().
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
int HICMA_zgytlr_Tile_Async(
        MORSE_enum       uplo,
        MORSE_desc_t     *AUV,
        MORSE_desc_t     *AD,
        MORSE_desc_t     *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int  compress_diag,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence,
        MORSE_request_t  *request)
{
    MORSE_desc_t *A = AUV; // FIXME

    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgytlr_Tile", "morse not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgytlr_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgytlr_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(A) != MORSE_SUCCESS) {
        morse_error("MORSE_zgytlr_Tile", "invalid descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb) {
        morse_error("MORSE_zgytlr_Tile", "only square tiles supported");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }*/

    /* Quick return */
    if (chameleon_min( A->m, A->n ) == 0)
        return MORSE_SUCCESS;

    hicma_pzgytlr(
            uplo,
            AUV,
            AD,
            Ark,
            seed,
            maxrank, tol,
            compress_diag,
            Dense,
            sequence, request);

    return MORSE_SUCCESS;
}
