/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/
/**
 * @file pzhagdm.c
 *
 *  This file contains parallel routines for generating dense matrices.
 *
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
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
#include <stdio.h>

/**
 * Generate a dense matrix.
 */
void hicma_pzhagdm(
        MORSE_enum uplo,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < Dense->mt; m++) {
        tempmm = m == Dense->mt-1 ? Dense->m-m*Dense->mb : Dense->mb;
        int ldam = BLKLDD(Dense, m);

        //for (n = 0; n < A->mt; n++) {
        //    tempnn = n == A->mt-1 ? A->m-n*A->mb : A->mb;
        for (n = 0; n < Dense->nt; n++) { //I hope this change does not break anything
            tempnn = n == Dense->nt-1 ? Dense->n-n*Dense->nb : Dense->nb;
            if(uplo == MorseLower && m < n)
                continue;
            else if(uplo == MorseUpper && m > n)
                continue;
            /*printf("%s %d: %d %d\n", __FILE__, __LINE__, m, n);*/
            HICMA_TASK_zhagdm(
                    &options,
                    tempmm, tempnn,
                    Dense,
                    ldam,
                    m, n
                    );
        }
    }
    RUNTIME_options_finalize(&options, morse);
    /*MORSE_TASK_dataflush_all(); removed in newer chameleon */
}
/**
 * Generate a dense matrix.
 * The diagonal tiles of problem are used.
 */
void hicma_pzhagdmdiag(
        MORSE_enum uplo,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < Dense->mt; m++) {
        tempmm = m == Dense->mt-1 ? Dense->m-m*Dense->mb : Dense->mb;
        int ldam = BLKLDD(Dense, m);

        /*printf("%s %d: %d %d\n", __FILE__, __LINE__, m, n);*/
        HICMA_TASK_zhagdmi(
                &options,
                tempmm, tempmm,
                Dense,
                ldam,
                m, 0,
                m, m //index of tile in problem
                );
    }
    RUNTIME_options_finalize(&options, morse);
    /*MORSE_TASK_dataflush_all(); removed in newer chameleon */
}
