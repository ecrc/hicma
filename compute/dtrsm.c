/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 *
 * @file dtrsm.c
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Kadir Akbudak
 * @date 2018-11-08
 *
 **/

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/

/**
 *
 * file dtrsm.c
 *
 * HICMA computational routines
 * HICMA is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for HICMA 1.0.0
 * author Jakub Kurzak
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * date 2010-11-15
 *
 **/

#include "hicma.h"
#include <stdio.h>
#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>


/***************************************************************************//**
 *
 *  HICMA_dtrsm_Tile - Computes triangular solve.
 *  Both A and B/X matrices are in Tile Low Rank (TLR) format.
 *  Tile equivalent of HICMA_dtrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = HicmaLeft:  A*X = B
 *          = HicmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = HicmaUpper: Upper triangle of A is stored;
 *          = HicmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = HicmaNoTrans:   A is transposed;
 *          = HicmaTrans:     A is not transposed;
 *          = HicmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = HicmaNonUnit: A is non unit;
 *          = HicmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = HicmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = HicmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = HicmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_dtrsm_Tile(HICMA_enum side, HICMA_enum uplo,
                      HICMA_enum transA, HICMA_enum diag,
                      double alpha, 
                      HICMA_desc_t *AUV,
                      HICMA_desc_t *AD,
                      HICMA_desc_t *Ark,
                      HICMA_desc_t *BUV,
                      HICMA_desc_t *Brk,
                      int rk,
                      int maxrk,
                      double acc
                      )
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s rk:%d maxrk:%d acc:%e alpha:%e\n",
                __LINE__, __func__,
                rk, maxrk, acc, alpha);
    }
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dtrsm_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);
    HICMA_dtrsm_Tile_Async(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, BUV, Brk, 
            rk, maxrk, acc,
            sequence, &request);
    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( AD, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    HICMA_Desc_Flush( BUV, sequence );
    HICMA_Desc_Flush( Brk, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(BUV);*/
    /*RUNTIME_desc_getoncpu(Brk); */
    
    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_dtrsm_Tile_Async - Computes triangular solve.
 *  Both A and B/X matrices are in Tile Low Rank (TLR) format.
 *  Non-blocking equivalent of HICMA_dtrsm_Tile().
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
 *******************************************************************************/
int HICMA_dtrsm_Tile_Async(HICMA_enum side, HICMA_enum uplo,
                            HICMA_enum transA, HICMA_enum diag,
                            double alpha, 
                            HICMA_desc_t *AUV,
                            HICMA_desc_t *AD,
                            HICMA_desc_t *Ark,
                            HICMA_desc_t *BUV,
                            HICMA_desc_t *Brk,
                            int rk,
                            int maxrk,
                            double acc,
                            HICMA_sequence_t *sequence, HICMA_request_t *request)
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s rk:%d maxrk:%d acc:%e alpha:%e\n",
                __LINE__, __func__,
                rk, maxrk, acc, alpha);
    }
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dtrsm_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dtrsm_Tile", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dtrsm_Tile", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (hicma_desc_check(AUV) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dtrsm_Tile", "invalid first descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    if (hicma_desc_check(BUV) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dtrsm_Tile", "invalid second descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb || B->nb != B->mb) {*/
        /*hicma_error("HiCMA_dtrsm_Tile", "only square tiles supported");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/
    if (side != HicmaLeft && side != HicmaRight) {
        hicma_error("HiCMA_dtrsm_Tile", "illegal value of side");
        return hicma_request_fail(sequence, request, -1);
    }
    if (uplo != HicmaUpper && uplo != HicmaLower) {
        hicma_error("HiCMA_dtrsm_Tile", "illegal value of uplo");
        return hicma_request_fail(sequence, request, -2);
    }
    if (transA != HicmaConjTrans && transA != HicmaNoTrans && transA != HicmaTrans) {
        hicma_error("HiCMA_dtrsm_Tile", "illegal value of transA");
        return hicma_request_fail(sequence, request, -3);
    }
    if (diag != HicmaUnit && diag != HicmaNonUnit) {
        hicma_error("HiCMA_dtrsm_Tile", "illegal value of diag");
        return hicma_request_fail(sequence, request, -4);
    }

    hicma_pdtrsm(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, BUV, Brk, 
            rk, maxrk, acc,
            sequence, request);

    return HICMA_SUCCESS;
}

/***************************************************************************//**
 *
 *  HICMA_dtrsm_Tile - Computes triangular solve.
 *  A matrix is in Tile Low Rank (TLR) format and B/X matrix is dense.
 *  Tile equivalent of HICMA_dtrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = HicmaLeft:  A*X = B
 *          = HicmaRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = HicmaUpper: Upper triangle of A is stored;
 *          = HicmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = HicmaNoTrans:   A is transposed;
 *          = HicmaTrans:     A is not transposed;
 *          = HicmaConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = HicmaNonUnit: A is non unit;
 *          = HicmaUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = HicmaUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = HicmaLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = HicmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_dtrsmd_Tile(HICMA_enum side, HICMA_enum uplo,
                      HICMA_enum transA, HICMA_enum diag,
                      double alpha, 
                      HICMA_desc_t *AUV,
                      HICMA_desc_t *AD,
                      HICMA_desc_t *Ark,
                      HICMA_desc_t *Bdense,
                      int maxrk
                      )
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s maxrk:%d alpha:%e\n",
                __LINE__, __func__,
                maxrk, alpha);
    }
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dtrsmd_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);
    HICMA_dtrsmd_Tile_Async(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, Bdense, 
            maxrk,
            sequence, &request);
    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( AD, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    HICMA_Desc_Flush( Bdense, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(Bdense);*/
    
    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_dtrsm_Tile_Async - Computes triangular solve.
 *  A matrix is in Tile Low Rank (TLR) format and B/X matrix is dense.
 *  Non-blocking equivalent of HICMA_dtrsm_Tile().
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
 *******************************************************************************/
int HICMA_dtrsmd_Tile_Async(HICMA_enum side, HICMA_enum uplo,
                            HICMA_enum transA, HICMA_enum diag,
                            double alpha, 
                            HICMA_desc_t *AUV,
                            HICMA_desc_t *AD,
                            HICMA_desc_t *Ark,
                            HICMA_desc_t *Bdense,
                            int maxrk,
                            HICMA_sequence_t *sequence, HICMA_request_t *request)
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s maxrk:%d alpha:%e\n",
                __LINE__, __func__,
                maxrk, alpha);
    }
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dtrsmd_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dtrsmd_Tile", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dtrsmd_Tile", "NULL request");
        return HICMA_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == HICMA_SUCCESS)
        request->status = HICMA_SUCCESS;
    else
        return hicma_request_fail(sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (hicma_desc_check(AUV) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dtrsmd_Tile", "invalid first descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    if (hicma_desc_check(Bdense) != HICMA_SUCCESS) {
        hicma_error("HiCMA_dtrsmd_Tile", "invalid second descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb || B->nb != B->mb) {*/
        /*hicma_error("HiCMA_dtrsm_Tile", "only square tiles supported");*/
        /*return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);*/
    /*}*/
    if (side != HicmaLeft && side != HicmaRight) {
        hicma_error("HiCMA_dtrsmd_Tile", "illegal value of side");
        return hicma_request_fail(sequence, request, -1);
    }
    if (uplo != HicmaUpper && uplo != HicmaLower) {
        hicma_error("HiCMA_dtrsmd_Tile", "illegal value of uplo");
        return hicma_request_fail(sequence, request, -2);
    }
    if (transA != HicmaConjTrans && transA != HicmaNoTrans && transA != HicmaTrans) {
        hicma_error("HiCMA_dtrsmd_Tile", "illegal value of transA");
        return hicma_request_fail(sequence, request, -3);
    }
    if (diag != HicmaUnit && diag != HicmaNonUnit) {
        hicma_error("HiCMA_dtrsmd_Tile", "illegal value of diag");
        return hicma_request_fail(sequence, request, -4);
    }

    hicma_pdtrsmd(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, Bdense, 
            maxrk,
            sequence, request);

    return HICMA_SUCCESS;
}

