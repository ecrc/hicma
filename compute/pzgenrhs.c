/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
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
 * author Rabab Alomairy
 * date 2020-02-15
 *
 **/
#include "morse.h"
#include "control/common.h"
#include "hicma_runtime_z.h"

/**
 *  hicma_pzgenrhs_virus - Generate a random matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  HICMA_zgenrhs_Tile - Generate RHS matrix by tiles.
 *
 *******************************************************************************
 *
 * @param[out] A
 *        Store RHS values in descriptor A  
 *
 * @param[in] sequence
 *          Identifies the sequence of function calls that this call belongs to
 *          (for completion checks and exception handling purposes).
 *
 * @param[out] request
 *          Identifies this function call (for exception handling purposes).
 *

 ******************************************************************************/

void hicma_pzgenrhs(
        MORSE_desc_t *A,
        MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n;
    int tempmm, tempnn;
    int index;
    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        int ldam = BLKLDD(A, m);
                if(1)HICMA_TASK_zgenrhs(
                        &options,
                        tempmm,
                        tempnn,
                        A, m, 0,
                        ldam,
                        A->m, m*A->mb, m*A->nb);
    }
    RUNTIME_options_finalize(&options, morse);
}
