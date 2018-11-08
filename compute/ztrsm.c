/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/
/**
 *
 * @file ztrsm.c
 *
 *  HICMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> s d c
 *
 **/
/*
 *
 * file ztrsm.c
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
 * precisions normal z -> s d c
 *
 */
#include "hicma.h"
#include <stdio.h>
#include "control/common.h"
#include "control/hicma_common.h"

/***************************************************************************//**
 *
 *  HICMA_ztrsm_Tile - Computes triangular solve.
 *  Both A and B/X matrices are in Tile Low Rank (TLR) format.
 *  Tile equivalent of HICMA_ztrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = MorseLeft:  A*X = B
 *          = MorseRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   A is transposed;
 *          = MorseTrans:     A is not transposed;
 *          = MorseConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = MorseNonUnit: A is non unit;
 *          = MorseUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = MorseUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = MorseLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = MorseUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_ztrsm_Tile(MORSE_enum side, MORSE_enum uplo,
                      MORSE_enum transA, MORSE_enum diag,
                      double alpha, 
                      MORSE_desc_t *AUV, 
                      MORSE_desc_t *AD, 
                      MORSE_desc_t *Ark, 
                      MORSE_desc_t *BUV,
                      MORSE_desc_t *Brk,
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
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsm_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    HICMA_ztrsm_Tile_Async(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, BUV, Brk, 
            rk, maxrk, acc,
            sequence, &request);
    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    MORSE_Desc_Flush( BUV, sequence );
    MORSE_Desc_Flush( Brk, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(BUV);*/
    /*RUNTIME_desc_getoncpu(Brk); */
    
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_ztrsm_Tile_Async - Computes triangular solve.
 *  Both A and B/X matrices are in Tile Low Rank (TLR) format.
 *  Non-blocking equivalent of HICMA_ztrsm_Tile().
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
int HICMA_ztrsm_Tile_Async(MORSE_enum side, MORSE_enum uplo,
                            MORSE_enum transA, MORSE_enum diag,
                            double alpha, 
                            MORSE_desc_t *AUV, 
                            MORSE_desc_t *AD, 
                            MORSE_desc_t *Ark, 
                            MORSE_desc_t *BUV,
                            MORSE_desc_t *Brk,
                            int rk,
                            int maxrk,
                            double acc,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s rk:%d maxrk:%d acc:%e alpha:%e\n",
                __LINE__, __func__,
                rk, maxrk, acc, alpha);
    }
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsm_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztrsm_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztrsm_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(AUV) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsm_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(BUV) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsm_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb || B->nb != B->mb) {*/
        /*morse_error("MORSE_ztrsm_Tile", "only square tiles supported");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/
    if (side != MorseLeft && side != MorseRight) {
        morse_error("MORSE_ztrsm_Tile", "illegal value of side");
        return morse_request_fail(sequence, request, -1);
    }
    if (uplo != MorseUpper && uplo != MorseLower) {
        morse_error("MORSE_ztrsm_Tile", "illegal value of uplo");
        return morse_request_fail(sequence, request, -2);
    }
    if (transA != MorseConjTrans && transA != MorseNoTrans && transA != MorseTrans) {
        morse_error("MORSE_ztrsm_Tile", "illegal value of transA");
        return morse_request_fail(sequence, request, -3);
    }
    if (diag != MorseUnit && diag != MorseNonUnit) {
        morse_error("MORSE_ztrsm_Tile", "illegal value of diag");
        return morse_request_fail(sequence, request, -4);
    }

    hicma_pztrsm(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, BUV, Brk, 
            rk, maxrk, acc,
            sequence, request);

    return MORSE_SUCCESS;
}

/***************************************************************************//**
 *
 *  HICMA_ztrsm_Tile - Computes triangular solve.
 *  A matrix is in Tile Low Rank (TLR) format and B/X matrix is dense.
 *  Tile equivalent of HICMA_ztrsm().
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          Specifies whether A appears on the left or on the right of X:
 *          = MorseLeft:  A*X = B
 *          = MorseRight: X*A = B
 *
 * @param[in] uplo
 *          Specifies whether the matrix A is upper triangular or lower triangular:
 *          = MorseUpper: Upper triangle of A is stored;
 *          = MorseLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or conjugate transposed:
 *          = MorseNoTrans:   A is transposed;
 *          = MorseTrans:     A is not transposed;
 *          = MorseConjTrans: A is conjugate transposed.
 *
 * @param[in] diag
 *          Specifies whether or not A is unit triangular:
 *          = MorseNonUnit: A is non unit;
 *          = MorseUnit:    A us unit.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = MorseUpper, the leading N-by-N upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = MorseLower, the leading N-by-N
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = MorseUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in,out] B
 *          On entry, the N-by-NRHS right hand side matrix B.
 *          On exit, if return value = 0, the N-by-NRHS solution matrix X.
 *
 *******************************************************************************
 *
 * @return
 *          \retval MORSE_SUCCESS successful exit
 *
 ******************************************************************************/
int HICMA_ztrsmd_Tile(MORSE_enum side, MORSE_enum uplo,
                      MORSE_enum transA, MORSE_enum diag,
                      double alpha, 
                      MORSE_desc_t *AUV, 
                      MORSE_desc_t *AD, 
                      MORSE_desc_t *Ark, 
                      MORSE_desc_t *Bdense,
                      int maxrk
                      )
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s maxrk:%d alpha:%e\n",
                __LINE__, __func__,
                maxrk, alpha);
    }
    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsmd_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);
    HICMA_ztrsmd_Tile_Async(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, Bdense, 
            maxrk,
            sequence, &request);
    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    MORSE_Desc_Flush( Bdense, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(Ark);*/
    /*RUNTIME_desc_getoncpu(Bdense);*/
    
    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/***************************************************************************//**
 *
 *  HICMA_ztrsm_Tile_Async - Computes triangular solve.
 *  A matrix is in Tile Low Rank (TLR) format and B/X matrix is dense.
 *  Non-blocking equivalent of HICMA_ztrsm_Tile().
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
int HICMA_ztrsmd_Tile_Async(MORSE_enum side, MORSE_enum uplo,
                            MORSE_enum transA, MORSE_enum diag,
                            double alpha, 
                            MORSE_desc_t *AUV, 
                            MORSE_desc_t *AD, 
                            MORSE_desc_t *Ark, 
                            MORSE_desc_t *Bdense,
                            int maxrk,
                            MORSE_sequence_t *sequence, MORSE_request_t *request)
{
    if(HICMA_get_print_index() == 1){
        printf("%d:%s maxrk:%d alpha:%e\n",
                __LINE__, __func__,
                maxrk, alpha);
    }
    MORSE_context_t *morse;

    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("MORSE_ztrsmd_Tile", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        morse_fatal_error("MORSE_ztrsmd_Tile", "NULL sequence");
        return MORSE_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        morse_fatal_error("MORSE_ztrsmd_Tile", "NULL request");
        return MORSE_ERR_UNALLOCATED;
    }
    /* Check sequence status */
    if (sequence->status == MORSE_SUCCESS)
        request->status = MORSE_SUCCESS;
    else
        return morse_request_fail(sequence, request, MORSE_ERR_SEQUENCE_FLUSHED);

    /* Check descriptors for correctness */
    if (morse_desc_check(AUV) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsmd_Tile", "invalid first descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    if (morse_desc_check(Bdense) != MORSE_SUCCESS) {
        morse_error("MORSE_ztrsmd_Tile", "invalid second descriptor");
        return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    /*if (A->nb != A->mb || B->nb != B->mb) {*/
        /*morse_error("MORSE_ztrsm_Tile", "only square tiles supported");*/
        /*return morse_request_fail(sequence, request, MORSE_ERR_ILLEGAL_VALUE);*/
    /*}*/
    if (side != MorseLeft && side != MorseRight) {
        morse_error("MORSE_ztrsmd_Tile", "illegal value of side");
        return morse_request_fail(sequence, request, -1);
    }
    if (uplo != MorseUpper && uplo != MorseLower) {
        morse_error("MORSE_ztrsmd_Tile", "illegal value of uplo");
        return morse_request_fail(sequence, request, -2);
    }
    if (transA != MorseConjTrans && transA != MorseNoTrans && transA != MorseTrans) {
        morse_error("MORSE_ztrsmd_Tile", "illegal value of transA");
        return morse_request_fail(sequence, request, -3);
    }
    if (diag != MorseUnit && diag != MorseNonUnit) {
        morse_error("MORSE_ztrsmd_Tile", "illegal value of diag");
        return morse_request_fail(sequence, request, -4);
    }

    hicma_pztrsmd(side, uplo, transA, diag, alpha, 
            AUV, AD, Ark, Bdense, 
            maxrk,
            sequence, request);

    return MORSE_SUCCESS;
}

