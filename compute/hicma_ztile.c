/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/

/**
 *
 * @file ztile.c
 *
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 */

#include <control/common.h>
#include <include/hicma_z.h>
#include <control/hicma_compute_z.h>
/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t
 *
 *  HICMA_zLapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 *******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[in,out] A
 *          Descriptor of the HICMA matrix in tile layout.
 *          If HICMA_TRANSLATION_MODE is set to HICMA_INPLACE,
 *          A->mat is not used and set to Af77 when returns, else if
 *          HICMA_TRANSLATION_MODE is set to HICMA_OUTOFPLACE,
 *          A->mat has to be allocated before.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa HICMA_zTile_to_Lapack
 * @sa HICMA_cLapack_to_Tile
 * @sa HICMA_dLapack_to_Tile
 * @sa HICMA_sLapack_to_Tile
 *
 */
int HICMA_zLapack_to_Tile( HICMA_Complex64_t *Af77, int LDA, HICMA_desc_t *A )
{
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request;
    HICMA_desc_t *B;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zLapack_to_Tile", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (hicma_desc_check( A ) != HICMA_SUCCESS) {
        hicma_error("HiCMA_zLapack_to_Tile", "invalid descriptor");
        return HICMA_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    HICMA_Desc_Create_User( &B, Af77, HicmaComplexDouble, A->mb, A->nb, A->bsiz,
                            LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                            hicma_getaddr_cm, hicma_getblkldd_cm, NULL );

    /* Start the computation */
    hicma_sequence_create( hicma, &sequence );

    hicma_pzlacpy( HicmaUpperLower, B, A, sequence, &request );

    HICMA_Desc_Flush( B, sequence );
    HICMA_Desc_Flush( A, sequence );

    hicma_sequence_wait( hicma, sequence );

    /* Destroy temporary B descriptor */
    HICMA_Desc_Destroy( &B );

    status = sequence->status;
    hicma_sequence_destroy( hicma, sequence );
    return status;
}

/**
 ********************************************************************************
 *
 * @ingroup HICMA_Complex64_t
 *
 *  HICMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 *******************************************************************************
 *
 * @param[in] A
 *          Descriptor of the HICMA matrix in tile layout.
 *
 * @param[in,out] Af77
 *          LAPACK matrix.
 *          If HICMA_TRANSLATION_MODE is set to HICMA_INPLACE,
 *          Af77 has to be A->mat, else if
 *          HICMA_TRANSLATION_MODE is set to HICMA_OUTOFPLACE,
 *          Af77 has to be allocated before.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 *******************************************************************************
 *
 * @sa HICMA_zLapack_to_Tile
 * @sa HICMA_cTile_to_Lapack
 * @sa HICMA_dTile_to_Lapack
 * @sa HICMA_sTile_to_Lapack
 *
 */
int HICMA_zTile_to_Lapack( HICMA_desc_t *A, HICMA_Complex64_t *Af77, int LDA )
{
    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request;
    HICMA_desc_t *B;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_zTile_to_Lapack", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    /* Check descriptor for correctness */
    if (hicma_desc_check( A ) != HICMA_SUCCESS) {
        hicma_error("HiCMA_zTile_to_Lapack", "invalid descriptor");
        return HICMA_ERR_ILLEGAL_VALUE;
    }

    /* Create the B descriptor to handle the Lapack format matrix */
    HICMA_Desc_Create_User( &B, Af77, HicmaComplexDouble, A->mb, A->nb, A->bsiz,
                            LDA, A->n, 0, 0, A->m, A->n, 1, 1,
                            hicma_getaddr_cm, hicma_getblkldd_cm, NULL );

    /* Start the computation */
    hicma_sequence_create( hicma, &sequence );

    hicma_pzlacpy( HicmaUpperLower, A, B, sequence, &request );

    HICMA_Desc_Flush( A, sequence );
    HICMA_Desc_Flush( B, sequence );

    hicma_sequence_wait( hicma, sequence );

    HICMA_Desc_Destroy( &B );

    status = sequence->status;
    hicma_sequence_destroy( hicma, sequence );
    return status;
}
