/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzgytlr.c
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
 * file pzgytlr.c
 *
 * MORSE auxiliary routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

#include <hicma.h>
#include <control/common.h>
#include <include/hicma_runtime_z.h>

extern int store_only_diagonal_tiles;

/***************************************************************************//**
 * Generate a spatial statistics matrix using STARS-H.
 **/
void hicma_pzgytlr(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        HICMA_desc_t *Dense,
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
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        int tempmmD = m == AD->mt-1 ? AD->m-m*AD->mb : AD->mb;
        int ldamUV = BLKLDD(AUV, m);
        int ldamD = BLKLDD(AD, m);

        //for (n = 0; n < A->mt; n++) {
        //    tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
        for (n = 0; n < A->nt; n++) { //I hope this change does not break anything
            if(0 && A->nt == 1) { //FIXME for B in TRSM
                tempnn = A->nb;
            } else {
                tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
            }

            // if(m<n)
            //     continue;
            if(uplo == HicmaLower && m < n)
                continue;
            else
            if(uplo == HicmaUpper && m > n)
                continue;
            //printf("Tile (%d,%d) ldam=%d A->m=%d A->n=%d A->mb=%d A->nb=%d\n ", m, n, ldam, A->m, A->n, A->mb, A->nb);

            int call_diag = 0;
            int AD_icol;
            if(store_only_diagonal_tiles == 1){
                if(m == n) {
                    call_diag = 1;
                    AD_icol = 0;
                } else {
                }
            } else {
                call_diag = 0; //FIXME call_diag must 1 when AD is full matrix and used
                AD_icol = n;
            }
            if(0 && AUV->nt == 1){
                call_diag = 0; //FIXME added for B in TRSM
            }
            if(call_diag == 1) {
                /*printf("diag %d,%d\n", m, n);*/
                HICMA_TASK_zgytlr_diag(
                        &options,
                        tempmmD, tempmmD,
                        AUV,
                        AD, m, AD_icol,
                        Ark,
                        m, n, 
                        ldamD, 0, 0,
                        A->m, m*A->mb, n*A->mb, seed,
                        maxrank, tol,
                        compress_diag,
                        Dense
                        );
            } else {
                //printf("off  %d,%d\n", m, n);
                HICMA_TASK_zgytlr(
                        &options,
                        tempmmD, tempnn,
                        AUV,
                        Ark,
                        m, n, 
                        ldamD, ldamUV, ldamUV,
                        A->m, m*A->mb, n*A->mb, seed,
                        maxrank, tol,
                        compress_diag,
                        Dense
                        );
            }
        }
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
    //HICMA_TASK_dataflush_all(); removed in newer chameleon


}
