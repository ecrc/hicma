/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 * @file dpotrf.c
 *
 * This file contains top-level functions for Cholesky factorization.
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/

/**
 * file dpotrf.c
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
 **/

#include <control/common.h>
#include <hicma_common.h>
#include <include/hicma_d.h>
#include <control/hicma_compute_d.h>


/***************************************************************************//**
 *
 *  HICMA_dpotrf_Tile - Computes the Cholesky factorization of a symmetric 
 *  positive definite matrix in tile low-rank (TLR) format.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *******************************************************************************
 *
 * @param[in] uplo
 *          = HicmaUpper: Upper triangle of A is stored (Not supported yet)
 *          = HicmaLower: Lower triangle of A is stored.
 *
 * @param[in] A
 *          On entry, the symmetric positive definite TLR matrix A.
 *          If uplo = HicmaUpper, the leading N-by-N upper triangular part of A
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
 *          \retval HICMA_SUCCESS successful exit
 *          \retval >0 if i, the leading minor of order i of A is not positive definite, so the
 *               factorization could not be completed, and the solution has not been computed.
 *
 ******************************************************************************/
int HICMA_dpotrf_Tile(HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        int rk, int maxrk, double acc
        )
{
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dpotrf_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);
    HICMA_dpotrf_Tile_Async(uplo,
            AUV, AD, Ark,
            rk, maxrk, acc,
            sequence, &request
            );
    HICMA_Desc_Flush( AD, sequence );
    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}
/***************************************************************************//**
 *
 *  HICMA_dpotrf_Tile_Async - Computes the Cholesky factorization of a symmetric
 *  positive definite positive definite matrix.
 *  Non-blocking equivalent of HICMA_dpotrf_Tile().
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
int HICMA_dpotrf_Tile_Async(HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        int rk, int maxrk, double acc,
        HICMA_sequence_t *sequence, HICMA_request_t *request
        )
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_dpotrf_Tile_Async", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_dpotrf_Tile_Async", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    if (request == NULL) {
        hicma_fatal_error("HiCMA_dpotrf_Tile_Async", "NULL request");
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
            || (hicma_desc_check(AD) != HICMA_SUCCESS)
            || (hicma_desc_check(Ark) != HICMA_SUCCESS)
            ){
        hicma_error("HiCMA_dpotrf_Tile_Async", "invalid descriptor");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    /* Check input arguments */
    if (AD->nb != AD->mb) {
        hicma_error("HiCMA_dpotrf_Tile_Async", "only square tiles supported");
        return hicma_request_fail(sequence, request, HICMA_ERR_ILLEGAL_VALUE);
    }
    if (uplo != HicmaUpper && uplo != HicmaLower) {
        hicma_error("HiCMA_dpotrf_Tile_Async", "illegal value of uplo");
        return hicma_request_fail(sequence, request, -1);
    }
    /* Quick return */
/*
    if (hicma_max(N, 0) == 0)
        return HICMA_SUCCESS;
*/

    hicma_pdpotrf(uplo, AUV, AD, Ark, sequence, request,
            rk, maxrk, acc
            );

    return HICMA_SUCCESS;
}
