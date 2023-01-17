/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzplrnt.c
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
 * file pzplrnt.c
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
 * @precisions normal z -> s d c
 *
 **/

#include <control/common.h>
#include <control/hicma_compute_z.h>
#include <include/hicma_runtime_z.h>

#define A(m, n) A,  m,  n

/**
 *  hicma_pzplghe - Generate a random matrix by tiles.
 */
void hicma_pzplrnt(HICMA_desc_t *A, unsigned long long int seed,
                         HICMA_sequence_t *sequence, HICMA_request_t *request) {
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int m, n;
    int ldam;
    int tempmm, tempnn;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
        ldam = BLKLDD(A, m);

        for (n = 0; n < A->nt; n++) {
            tempnn = n == A->nt - 1 ? A->n - n * A->nb : A->nb;

            HICMA_TASK_zplrnt(
                    &options,
                    tempmm, tempnn, A(m, n), ldam,
                    A->m, m * A->mb, n * A->nb, seed);
        }
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
}
