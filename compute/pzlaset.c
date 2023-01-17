/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzlaset.c
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
 * file pzlaset.c
 *
 * @brief Chameleon zlaset parallel algorithm
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for HICMA 1.0.0
 * @author Hatem Ltaief
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

/**
 *  Parallel initialization a 2-D array A to BETA on the diagonal and
 *  ALPHA on the offdiagonals.
 */
void hicma_pzlaset(HICMA_enum uplo,
                         HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                         HICMA_desc_t *A,
                         HICMA_sequence_t *sequence, HICMA_request_t *request) {
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int i, j;
    int ldai, ldaj;
    int tempim;
    int tempjm, tempjn;
    int minmn = hicma_min(A->mt, A->nt);

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;

    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    if (uplo == HicmaLower) {
        for (j = 0; j < minmn; j++) {
            tempjm = j == A->mt - 1 ? A->m - j * A->mb : A->mb;
            tempjn = j == A->nt - 1 ? A->n - j * A->nb : A->nb;
            ldaj = BLKLDD(A, j);
            HICMA_TASK_zlaset(
                    &options,
                    HicmaLower, tempjm, tempjn, alpha, beta,
                    A(j, j), ldaj);

            for (i = j + 1; i < A->mt; i++) {
                tempim = i == A->mt - 1 ? A->m - i * A->mb : A->mb;
                ldai = BLKLDD(A, i);
                HICMA_TASK_zlaset(
                        &options,
                        HicmaUpperLower, tempim, tempjn, alpha, alpha,
                        A(i, j), ldai);
            }
        }
    } else if (uplo == HicmaUpper) {
        for (j = 1; j < A->nt; j++) {
            tempjn = j == A->nt - 1 ? A->n - j * A->nb : A->nb;
            for (i = 0; i < hicma_min(j, A->mt); i++) {
                tempim = i == A->mt - 1 ? A->m - i * A->mb : A->mb;
                ldai = BLKLDD(A, i);
                HICMA_TASK_zlaset(
                        &options,
                        HicmaUpperLower, tempim, tempjn, alpha, alpha,
                        A(i, j), ldai);
            }
        }
        for (j = 0; j < minmn; j++) {
            tempjm = j == A->mt - 1 ? A->m - j * A->mb : A->mb;
            tempjn = j == A->nt - 1 ? A->n - j * A->nb : A->nb;
            ldaj = BLKLDD(A, j);
            HICMA_TASK_zlaset(
                    &options,
                    HicmaUpper, tempjm, tempjn, alpha, beta,
                    A(j, j), ldaj);
        }
    } else {
        for (i = 0; i < A->mt; i++) {
            tempim = i == A->mt - 1 ? A->m - i * A->mb : A->mb;
            ldai = BLKLDD(A, i);
            for (j = 0; j < A->nt; j++) {
                tempjn = j == A->nt - 1 ? A->n - j * A->nb : A->nb;
                HICMA_TASK_zlaset(
                        &options,
                        HicmaUpperLower, tempim, tempjn, alpha, alpha,
                        A(i, j), ldai);
            }
        }
        for (j = 0; j < minmn; j++) {
            tempjm = j == A->mt - 1 ? A->m - j * A->mb : A->mb;
            tempjn = j == A->nt - 1 ? A->n - j * A->nb : A->nb;
            ldaj = BLKLDD(A, j);
            HICMA_TASK_zlaset(
                    &options,
                    HicmaUpperLower, tempjm, tempjn, alpha, beta,
                    A(j, j), ldaj);
        }
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
}
