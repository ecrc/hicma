/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
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
 * @file pdgenrhs.c
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @author Rabab Alomairy
 * @date 2020-02-15
 *
 **/

#include <hicma.h>
#include <control/common.h>
#include <include/hicma_runtime_d.h>
#include <control/hicma_compute_d.h>

/**
 *  hicma_pdgenrhs_virus - Generate a random matrix by tiles.
 *  Operates on matrices stored by tiles.
 *  All matrices are passed through descriptors.
 *  All dimensions are taken from the descriptors.
 *
 *  HICMA_dgenrhs_Tile - Generate RHS matrix by tiles.
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

void hicma_pdgenrhs(
        HICMA_desc_t *A,
        HICMA_sequence_t *sequence, HICMA_request_t *request )
{
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int m, n;
    int tempmm, tempnn;
    int index;
    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    for (m = 0; m < A->mt; m++) {
        tempmm = m == A->mt-1 ? A->m-m*A->mb : A->mb;
        tempnn = n == A->nt-1 ? A->n-n*A->nb : A->nb;
        int ldam = BLKLDD(A, m);
                HICMA_TASK_dgenrhs(
                        &options,
                        tempmm,
                        tempnn,
                        A, m, 0,
                        ldam,
                        A->m, m*A->mb, m*A->nb);
    }
    HICMA_RUNTIME_options_finalize(&options, hicma);
}
