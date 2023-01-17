/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 * @file ddiag.c
 *
 * This file contains the function for copying tiles of a tile vector into diagonal tiles of a tile matrix.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2019-11-14
 **/

#include <hicma.h>
#include <assert.h>

#include <control/common.h>
#include <include/hicma_runtime_d.h>

extern int store_only_diagonal_tiles;

/*
 * Uncompresses lower triangular part. Computes D=U*V^T. Ranks of U and Vs stored in Ark
 */
int HICMA_ddiag_vec2mat(
        HICMA_desc_t *vec, HICMA_desc_t *mat) {

    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;
    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_diag_vec2mat", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);


    /*HICMA_context_t *hicma;*/
    HICMA_option_t options;
    /*hicma = hicma_context_self();*/
    if (sequence->status != HICMA_SUCCESS)
        return HICMA_ERR_NOT_INITIALIZED;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, &request);
    assert(vec->mb == mat->mb);
    assert(vec->nb == mat->nb);
    assert(vec->mb == vec->nb);
    int i;
    for (i = 0; i < vec->mt; i++) {
        int vecicol;
        if (store_only_diagonal_tiles == 1) {
            vecicol = 0;
        } else {
            vecicol = i;
        }
        //@KADIR FIXME handle leftovers
        int ldv = BLKLDD(vec, i);
        int ldm = BLKLDD(mat, i);
        int tempii = i == vec->mt - 1 ? vec->m - i * vec->mb : vec->mb;
        //printf("i=%d ldv=%d ldm=%d vec->mb=%d mat->mb=%d tempii=%d\n", i, ldv, ldm, vec->mb, mat->mb, tempii);
        HICMA_TASK_dlacpy( //FIXME convert to z
                &options,
                HicmaUpperLower,
                tempii, tempii, vec->mb,
                vec, i, vecicol, ldv,
                mat, i, i, ldm);
    }
    HICMA_RUNTIME_sequence_wait(hicma, sequence);
    HICMA_RUNTIME_options_finalize(&options, hicma);
    //HICMA_TASK_dataflush_all(); removed in newer chameleon

    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    //HICMA_RUNTIME_options_finalize(&options, hicma);


    HICMA_Desc_Flush(vec, sequence);
    HICMA_Desc_Flush(mat, sequence);

    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(vec);*/
    /*RUNTIME_desc_getoncpu(mat);*/

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}
