/*
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/


/**
 * @file zgemm.c
 *
 * This file contains top-level functions for matrix-matrix multiplication.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
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

/*
 *
 * @file zgemm.c
 *
 *  MORSE computational routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2018-11-08
 * @precisions normal z -> s d c
 *
 **/
#include "morse.h"
#include "control/common.h"
#include "control/hicma_common.h"


/***************************************************************************//**
 *
 * HICMA_zgemm_Tile - Performs  multiplication of tile-low-rank (TLR) matrices
 * in the form of 
 *    \f[ C = \alpha [op( A )\times op( B )] + \beta C \f],
 * where op( X ) is one of
 *
 *    op( X ) = X  or op( X ) = X' or op( X ) = conjg( X' )
 *
 * alpha and beta are scalars, and A, B and C which are all in TLR format.
 * Operates on matrices stored by tiles.
 * All matrices are passed through descriptors.
 * All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed or not transposed:
 *          = MorseNoTrans:   A is not transposed;
 *          = MorseTrans:     A is transposed;
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed or not transposed:
 *          = MorseNoTrans:   B is not transposed;
 *          = MorseTrans:     B is transposed;
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] AUV
 *          AUV is a descriptor for the tile low-rank matrix A 
 *
 * @param[in] Ark
 *          Ark is a descriptor containing rank of each tile of A
 *
 * @param[in] BUV
 *          BUV is a descriptor for the tile low-rank matrix B 
 *
 * @param[in] Brk
 *          Brk is a descriptor containing rank of each tile of B
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[out] CUV
 *          CUV is a descriptor for the tile low-rank matrix C.
 *          On exit, it is over-written by the result TLR matrix.
 *
 * @param[out] Crk
 *          Crk is a descriptor containing rank of each tile of C.
 *          On exit, it is over-written by the new rank of each tile. 
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_zgemm_Tile(MORSE_enum transA, MORSE_enum transB,
        double alpha,
        MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        MORSE_desc_t *BUV, MORSE_desc_t *Brk,
        double beta,
        MORSE_desc_t *CUV, MORSE_desc_t *Crk,
        int rk,
        int maxrk,
        double acc
        )
{
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);

    HICMA_zgemm_Tile_Async( transA, transB,
                            alpha, AUV, Ark,
                                   BUV, Brk,
                            beta,  CUV, Crk,
                            rk, maxrk, acc,
                            sequence, &request
                            );

    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( BUV, sequence );
    MORSE_Desc_Flush( CUV, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    MORSE_Desc_Flush( Brk, sequence );
    MORSE_Desc_Flush( Crk, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(BUV);*/
    /*RUNTIME_desc_getoncpu(CUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(Brk);*/
    /*RUNTIME_desc_getoncpu(Crk);*/

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_zgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of HICMA_zgemm_Tile().
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
int HICMA_zgemm_Tile_Async(MORSE_enum transA, MORSE_enum transB,
        double alpha,
        MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        MORSE_desc_t *BUV, MORSE_desc_t *Brk,
        double beta,
        MORSE_desc_t *CUV, MORSE_desc_t *Crk,
        int rk,
        int maxrk,
        double acc ,
        MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    MORSE_context_t *morse;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_zgemm_Tile_Async", "NULL request");
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
        ||(morse_desc_check(BUV) != MORSE_SUCCESS)
        ||(morse_desc_check(CUV) != MORSE_SUCCESS)
        ||(morse_desc_check(Ark) != MORSE_SUCCESS)
        ||(morse_desc_check(Brk) != MORSE_SUCCESS)
        ||(morse_desc_check(Crk) != MORSE_SUCCESS)
        )
    {
        morse_error("MORSE_zgemm_Tile_Async", "some invalid descriptors");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((transA != MorseNoTrans) && (transA != MorseTrans) && (transA != MorseConjTrans)) {
        morse_error("MORSE_zgemm_Tile_Async", "illegal value of transA");
        return morse_request_fail(sequence, request, -1);
    }
    if ((transB != MorseNoTrans) && (transB != MorseTrans) && (transB != MorseConjTrans)) {
        morse_error("MORSE_zgemm_Tile_Async", "illegal value of transB");
        return morse_request_fail(sequence, request, -2);
    }

    if ( transA == MorseNoTrans ) {
        Am  = AUV->m;
        An  = AUV->n;
        Amb = AUV->mb;
        Anb = AUV->nb;
        Ai  = AUV->i;
        Aj  = AUV->j;
    } else {
        Am  = AUV->n;
        An  = AUV->m;
        Amb = AUV->nb;
        Anb = AUV->mb;
        Ai  = AUV->j;
        Aj  = AUV->i;
    }

    if ( transB == MorseNoTrans ) {
        Bm  = BUV->m;
        Bn  = BUV->n;
        Bmb = BUV->mb;
        Bnb = BUV->nb;
        Bi  = BUV->i;
        Bj  = BUV->j;
    } else {
        Bm  = BUV->n;
        Bn  = BUV->m;
        Bmb = BUV->nb;
        Bnb = BUV->mb;
        Bi  = BUV->j;
        Bj  = BUV->i;
    }

    // Commented out by @kadir because Currently, Hicma's tiles can be tall and skinny matrices.
    // these conditions need more validations
    /*if ( (Amb != CUV->mb) || (Anb != Bmb) || (Bnb != CUV->nb) ) {*/
        /*morse_error("MORSE_zgemm_Tile_Async", "tile sizes have to match");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ( (Am != CUV->m) || (An != Bm) || (Bn != CUV->n) ) {*/
        /*morse_error("MORSE_zgemm_Tile_Async", "sizes of matrices have to match");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ( (Ai != CUV->i) || (Aj != Bi) || (Bj != CUV->j) ) {*/
        /*morse_error("MORSE_zgemm_Tile_Async", "start indexes have to match");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ((AUV->nb != AUV->mb) || (BUV->nb != BUV->mb) || (CUV->nb != CUV->mb) ){*/
        /*morse_error("HICMA_zpotrf_Tile_Async", "only square tiles supported");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/

    M = CUV->m;
    N = CUV->n;
    K = An;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == 0.0 || K == 0) && beta == 1.0))
        // ((alpha == (MORSE_Complex64_t)0.0 || K == 0) && beta == (MORSE_Complex64_t)1.0))
        return MORSE_SUCCESS;

    hicma_pzgemm(transA, transB,
                 alpha, AUV, Ark,
                        BUV, Brk,
                 beta,  CUV, Crk,
                 sequence, request,
                 rk, maxrk, acc);

    return MORSE_SUCCESS;
}
