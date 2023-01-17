/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pdhagcm.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
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
 */

/**
 *
 * file pdhagcm.c
 *
 * MORSE auxiliary routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

#include <hicma.h>
#include <control/common.h>
#include "hicma_runtime_d.h"
extern int store_only_diagonal_tiles;

/**
 * Generate a compressed matrix.
 * HicmaLower and HicmaUpper do not include diagnal tiles.
 **/
void hicma_pdhagcm(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank, double tol,
        HICMA_sequence_t *sequence, HICMA_request_t *request )
{
    HICMA_desc_t *A = AUV; // FIXME
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int m, n;
    int tempmm, tempnn;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? numrows_matrix-m*numrows_block : numrows_block;
        int ldamUV = BLKLDD(AUV, m);

        //for (n = 0; n < A->mt; n++) {
        //    tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
        for (n = 0; n < A->nt; n++) { //I hope this change does not break anything
            tempnn = n == A->nt-1 ? numcols_matrix-n*numcols_block : numcols_block;

            // if(m<n)
            //     continue;
            if(uplo == HicmaLower && m <= n)
                continue;
            else
            if(uplo == HicmaUpper && m >= n)
                continue;
            //printf("Tile (%d,%d) ldamUV=%d A->m=%d A->n=%d A->mb=%d A->nb=%d tempmm=%d tempnn=%d (%dx%d) (%dx%d) (%dx%d)\n", m, n, ldamUV, A->m, A->n, A->mb, A->nb, tempmm, tempnn, numrows_matrix, numcols_matrix, numrows_block, numcols_block, A->mt, A->nt);

            HICMA_TASK_dhagcm(
                    &options,
                    tempmm, tempnn,
                    AUV,
                    Ark,
                    m, n, 
                    ldamUV, ldamUV,
                    maxrank, tol, A->mt
                    );
        }
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
    //HICMA_TASK_dataflush_all(); removed in newer chameleon
}
