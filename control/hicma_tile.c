/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 *
 * @file tile.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon layout conversion wrappers
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 ***
 *
 * @defgroup Tile
 * @brief Group routines exposed to users for matrices conversion LAPACK-Tile
 *
 */

#include <control/hicma_auxiliary.h>
#include <hicma_common.h>
#include <include/hicma_z.h>
#include <include/hicma_d.h>
#include <include/hicma_c.h>
#include <include/hicma_s.h>
/**
 *
 * @ingroup Tile
 *
 *  HICMA_Lapack_to_Tile - Conversion from LAPACK layout to tile layout.
 *
 ******************************************************************************
 *
 * @param[in] Af77
 *          LAPACK matrix.
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 * @param[out] A
 *          Descriptor of the HICMA matrix in tile layout.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Lapack_to_Tile(void *Af77, int LDA, HICMA_desc_t *A)
{
    switch( A->dtyp ) {
        case HicmaComplexDouble:
            return HICMA_zLapack_to_Tile( (HICMA_Complex64_t *)Af77, LDA, A );
            break;
        case HicmaComplexFloat:
            return HICMA_cLapack_to_Tile( (HICMA_Complex32_t *)Af77, LDA, A );
            break;
        case HicmaRealFloat:
            return HICMA_sLapack_to_Tile( (float *)Af77, LDA, A );
            break;
        case HicmaRealDouble:
        default:
            return HICMA_dLapack_to_Tile( (double *)Af77, LDA, A );
    }
    return HICMA_ERR_ILLEGAL_VALUE;
}

/**
 *
 * @ingroup Tile
 *
 *  HICMA_Tile_to_Lapack - Conversion from tile layout to LAPACK layout.
 *
 ******************************************************************************
 *
 * @param[out] A
 *          Descriptor of the HICMA matrix in tile layout.
 *
 * @param[in] Af77
 *          LAPACK matrix (only needed on proc 0).
 *
 * @param[in] LDA
 *          The leading dimension of the matrix Af77.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Tile_to_Lapack(HICMA_desc_t *A, void *Af77, int LDA)
{
    switch( A->dtyp ) {
        case HicmaComplexDouble:
            return HICMA_zTile_to_Lapack( A, (HICMA_Complex64_t *)Af77, LDA );
            break;
        case HicmaComplexFloat:
            return HICMA_cTile_to_Lapack( A, (HICMA_Complex32_t *)Af77, LDA );
            break;
        case HicmaRealFloat:
            return HICMA_sTile_to_Lapack( A, (float *)Af77, LDA );
            break;
        case HicmaRealDouble:
        default:
            return HICMA_dTile_to_Lapack( A, (double *)Af77, LDA );
    }
    return HICMA_ERR_ILLEGAL_VALUE;
}
