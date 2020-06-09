/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzgytlr.c
 *
 * HiCMA auxiliary routines
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
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
 *
 * file pzplrnt.c
 *
 * MORSE auxiliary routines
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
 *
 **/
#include "morse.h"
#include "control/common.h"
#include "hicma_runtime_z.h"
extern int store_only_diagonal_tiles;

/***************************************************************************//**
 * Generate a spatial statistics matrix using STARS-H.
 **/
void hicma_pzcpy(
        MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AUVnew,
        MORSE_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_desc_t *A = AUV; // FIXME
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int tempmm, tempnn, tempnnnew;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        int tempmmAUVnew = m == AUVnew->mt-1 ? AUVnew->m-m*AUVnew->mb : AUVnew->mb;
        int ldamUV = BLKLDD(AUV, m);
        int ldamUVnew = BLKLDD(AUVnew, m);

        //for (n = 0; n < A->mt; n++) {
        //    tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
        for (n = 0; n < A->nt; n++) { //I hope this change does not break anything
                tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
                tempnnnew = n == A->nt-1 ? A->n-m*A->nb : A->nb;
            // if(m<n)
            //     continue;
            if(uplo == MorseLower && m < n)
                continue;
            else
            if(uplo == MorseUpper && m > n)
                continue;
            //printf("Tile (%d,%d) ldam=%d A->m=%d A->n=%d A->mb=%d A->nb=%d\n ", m, n, ldam, A->m, A->n, A->mb, A->nb);

            if(m!=n){
                //printf("off  %d,%d\n", m, n);
                HICMA_TASK_zcpy(
                        &options,
                        tempmm, tempnn,
                        AUV,
                        AUVnew,
                        Ark,
                        m, n, 
                        ldamUV, ldamUVnew, 
                        A->m, m*A->mb, n*A->mb, seed,
                        maxrank, tol,
                        compress_diag
                        );
            }
        }
    }
    RUNTIME_options_finalize(&options, morse);
    //MORSE_TASK_dataflush_all(); removed in newer chameleon


}
