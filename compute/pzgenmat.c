/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzgenmat.c 
 *
 * HiCMA auxiliary routines
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Rabab Alomairy
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
 * Generate a application matrix using STARS-H.
 **/
void hicma_pzgenmat(
        MORSE_desc_t *A,
        MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int i, j;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (i = 0; i < A->mt; i++) {
        int lda = BLKLDD(A, i);
        tempmm = i == A->mt-1 ? A->m-i*A->mb : A->mb;
        for (j = 0; j < A->nt; j++) {
            tempnn = j == A->nt-1 ? A->n-j*A->nb : A->nb;
            HICMA_TASK_zgenmat(&options, A, lda, i, j, tempmm, tempnn); 
            }
        }
    RUNTIME_options_finalize(&options, morse);
}
