/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file zdiag.c
 *
 * This file contains the function for copying tiles of a tile vector into diagonal tiles of a tile matrix.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#include "morse.h"
#include "control/common.h"
#include "hicma_runtime_z.h"
#include <assert.h>
extern int store_only_diagonal_tiles;
/*
 * Uncompresses lower triangular part. Computes D=U*V^T. Ranks of U and Vs stored in Ark
 */
int HICMA_zdiag_vec2mat(
        MORSE_desc_t *vec, MORSE_desc_t *mat)
{

    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;
    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_diag_vec2mat", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);


    /*MORSE_context_t *morse;*/
    MORSE_option_t options;
    /*morse = morse_context_self();*/
    if (sequence->status != MORSE_SUCCESS)
        return MORSE_ERR_NOT_INITIALIZED;
    RUNTIME_options_init(&options, morse, sequence, &request);
    assert(vec->mb == mat->mb);
    assert(vec->nb == mat->nb);
    assert(vec->mb == vec->nb);
    double zzero = (double) 0.0;
    double zone  = (double) 1.0;
    int i;
    for (i = 0; i < vec->mt; i++) {
        int vecicol;
        if(store_only_diagonal_tiles == 1){
            vecicol = 0;
        } else {
            vecicol = i;
        }
        //@KADIR FIXME handle leftovers
        int ldv = BLKLDD(vec, i);
        int ldm = BLKLDD(mat, i);
        int tempii = i == vec->mt-1 ? vec->m-i*vec->mb : vec->mb;
        //printf("i=%d ldv=%d ldm=%d vec->mb=%d mat->mb=%d tempii=%d\n", i, ldv, ldm, vec->mb, mat->mb, tempii);
        MORSE_TASK_dlacpy( //FIXME convert to z
                &options,
                MorseUpperLower,
                tempii, tempii, vec->mb,
                vec, i, vecicol, ldv,
                mat, i, i, ldm );
    }
    RUNTIME_sequence_wait( morse, sequence );
    RUNTIME_options_finalize( &options, morse );
    //MORSE_TASK_dataflush_all(); removed in newer chameleon

    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    //RUNTIME_options_finalize(&options, morse);


    MORSE_Desc_Flush( vec, sequence );
    MORSE_Desc_Flush( mat, sequence );

    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(vec);*/
    /*RUNTIME_desc_getoncpu(mat);*/

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}
