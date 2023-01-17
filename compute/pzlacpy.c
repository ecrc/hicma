/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzlacpy.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
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
 * file pzlacpy.c
 *
 * @brief Chameleon zlacpy parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> s d c
 *
 **/

#include <control/common.h>
#include <include/hicma_runtime_z.h>

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n
void hicma_pzlacpy(HICMA_enum uplo, HICMA_desc_t *A, HICMA_desc_t *B,
                         HICMA_sequence_t *sequence, HICMA_request_t *request) {
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int X, Y;
    int m, n;
    int ldam, ldbm;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    switch (uplo) {
        /*
         *  HicmaUpper
         */
        case HicmaUpper:
            for (m = 0; m < A->mt; m++) {
                X = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);
                if (m < A->nt) {
                    Y = m == A->nt - 1 ? A->n - m * A->nb : A->nb;
                    HICMA_TASK_zlacpy(
                            &options,
                            HicmaUpper,
                            X, Y, A->mb,
                            A(m, m), ldam,
                            B(m, m), ldbm);
                }
                for (n = m + 1; n < A->nt; n++) {
                    Y = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
                    HICMA_TASK_zlacpy(
                            &options,
                            HicmaUpperLower,
                            X, Y, A->mb,
                            A(m, n), ldam,
                            B(m, n), ldbm);
                }
            }
            break;
            /*
             *  HicmaLower
             */
        case HicmaLower:
            for (m = 0; m < A->mt; m++) {
                X = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);
                if (m < A->nt) {
                    Y = m == A->nt - 1 ? A->n - m * A->nb : A->nb;
                    HICMA_TASK_zlacpy(
                            &options,
                            HicmaLower,
                            X, Y, A->mb,
                            A(m, m), ldam,
                            B(m, m), ldbm);
                }
                for (n = 0; n < hicma_min(m, A->nt); n++) {
                    Y = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
                    HICMA_TASK_zlacpy(
                            &options,
                            HicmaUpperLower,
                            X, Y, A->mb,
                            A(m, n), ldam,
                            B(m, n), ldbm);
                }
            }
            break;
            /*
             *  HicmaUpperLower
             */
        case HicmaUpperLower:
        default:
            for (m = 0; m < A->mt; m++) {
                X = m == A->mt - 1 ? A->m - m * A->mb : A->mb;
                ldam = BLKLDD(A, m);
                ldbm = BLKLDD(B, m);
                for (n = 0; n < A->nt; n++) {
                    Y = n == A->nt - 1 ? A->n - n * A->nb : A->nb;
                    HICMA_TASK_zlacpy(
                            &options,
                            HicmaUpperLower,
                            X, Y, A->mb,
                            A(m, n), ldam,
                            B(m, n), ldbm);
                }
            }
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
}
