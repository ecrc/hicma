/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file dgemm.c
 *
 * This file contains top-level functions for matrix-matrix multiplication.
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
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
 * @file dgemm.c
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
 *
 **/

#include <hicma.h>
#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>


/***************************************************************************//**
 *
 * HICMA_dgemm_Tile - Performs  multiplication of tile-low-rank (TLR) matrices
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
 *          = HicmaNoTrans:   A is not transposed;
 *          = HicmaTrans:     A is transposed;
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed or not transposed:
 *          = HicmaNoTrans:   B is not transposed;
 *          = HicmaTrans:     B is transposed;
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
 *          \retval HICMA_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_dgemm_Tile(HICMA_enum transA, HICMA_enum transB,
        double alpha,
        HICMA_desc_t *AUV, HICMA_desc_t *Ark,
        HICMA_desc_t *BUV, HICMA_desc_t *Brk,
        double beta,
        HICMA_desc_t *CUV, HICMA_desc_t *Crk,
        int rk,
        int maxrk,
        double acc
        )
{
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dgemm_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);

    HICMA_dgemm_Tile_Async( transA, transB,
                            alpha, AUV, Ark,
                                   BUV, Brk,
                            beta,  CUV, Crk,
                            rk, maxrk, acc,
                            sequence, &request
                            );

    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( BUV, sequence );
    HICMA_Desc_Flush( CUV, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    HICMA_Desc_Flush( Brk, sequence );
    HICMA_Desc_Flush( Crk, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(BUV);*/
    /*RUNTIME_desc_getoncpu(CUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(Brk);*/
    /*RUNTIME_desc_getoncpu(Crk);*/

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_dgemm_Tile_Async - Performs matrix multiplication.
 *  Non-blocking equivalent of HICMA_dgemm_Tile().
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
int HICMA_dgemm_Tile_Async(HICMA_enum transA, HICMA_enum transB,
        double alpha,
        HICMA_desc_t *AUV, HICMA_desc_t *Ark,
        HICMA_desc_t *BUV, HICMA_desc_t *Brk,
        double beta,
        HICMA_desc_t *CUV, HICMA_desc_t *Crk,
        int rk,
        int maxrk,
        double acc ,
        HICMA_sequence_t *sequence, HICMA_request_t *request)
{
    HICMA_context_t *hicma;
    int M, N, K;
    int Am, An, Ai, Aj, Amb, Anb;
    int Bm, Bn, Bi, Bj, Bmb, Bnb;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dgemm_Tile_Async", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dgemm_Tile_Async", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dgemm_Tile_Async", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (
          (hicma_desc_check(AUV) != HICMA_SUCCESS)
        ||(hicma_desc_check(BUV) != HICMA_SUCCESS)
        ||(hicma_desc_check(CUV) != HICMA_SUCCESS)
        ||(hicma_desc_check(Ark) != HICMA_SUCCESS)
        ||(hicma_desc_check(Brk) != HICMA_SUCCESS)
        ||(hicma_desc_check(Crk) != HICMA_SUCCESS)
        )
    {
        hicma_error("HiCMA_dgemm_Tile_Async", "some invalid descriptors");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if ((transA != HicmaNoTrans) && (transA != HicmaTrans) && (transA != HicmaConjTrans)) {
        hicma_error("HiCMA_dgemm_Tile_Async", "illegal value of transA");
        return hicma_request_fail(sequence, request, -1);
    }
    if ((transB != HicmaNoTrans) && (transB != HicmaTrans) && (transB != HicmaConjTrans)) {
        hicma_error("HiCMA_dgemm_Tile_Async", "illegal value of transB");
        return hicma_request_fail(sequence, request, -2);
    }

    if ( transA == HicmaNoTrans ) {
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

    if ( transB == HicmaNoTrans ) {
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
        /*hicma_error("HiCMA_dgemm_Tile_Async", "tile sizes have to match");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ( (Am != CUV->m) || (An != Bm) || (Bn != CUV->n) ) {*/
        /*hicma_error("HiCMA_dgemm_Tile_Async", "sizes of matrices have to match");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ( (Ai != CUV->i) || (Aj != Bi) || (Bj != CUV->j) ) {*/
        /*hicma_error("HiCMA_dgemm_Tile_Async", "start indexes have to match");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/
    /*if ((AUV->nb != AUV->mb) || (BUV->nb != BUV->mb) || (CUV->nb != CUV->mb) ){*/
        /*hicma_error("HiCMA_dpotrf_Tile_Async", "only square tiles supported");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/

    M = CUV->m;
    N = CUV->n;
    K = An;

    /* Quick return */
    if (M == 0 || N == 0 ||
        ((alpha == 0.0 || K == 0) && beta == 1.0))
        // ((alpha == (HICMA_Complex64_t)0.0 || K == 0) && beta == (HICMA_Complex64_t)1.0))
        return HICMA_SUCCESS;

    hicma_pdgemm(transA, transB,
                 alpha, AUV, Ark,
                        BUV, Brk,
                 beta,  CUV, Crk,
                 sequence, request,
                 rk, maxrk, acc);

    return HICMA_SUCCESS;
}
