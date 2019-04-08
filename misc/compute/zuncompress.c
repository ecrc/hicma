/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * This file contains the function for uncompressing a tile low-rank matrix. 
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
/**
 * Uncompresses tile low-rank matrix AUV in to AD. 
 * Computes D=U*V^T. Ranks of U and Vs stored in Ark
 */
int HICMA_zuncompress(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark)
{

    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;
    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_zuncompress", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);


    /*MORSE_context_t *morse;*/
    MORSE_option_t options;
    /*morse = morse_context_self();*/
    if (sequence->status != MORSE_SUCCESS)
        return MORSE_ERR_NOT_INITIALIZED;
    RUNTIME_options_init(&options, morse, sequence, &request);
    double zzero = (double) 0.0;
    double zone  = (double) 1.0;
    int i, j;
    for (i = 0; i < AD->mt; i++) {
        /**
         * tempi contains number of rows of tiles in ith row. 
         */
        int tempi = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
        for (j = 0; j < AD->mt; j++) {
            if(uplo == MorseLower && i<=j)
                continue;
            else if(uplo == MorseUpper && i>=j)
                continue;
            //printf("%s %d: i:%d j:%d\n", __FILE__,__LINE__,i, j);
            /**
             * leading dimension of AUV. It might be equal to tempi.
             * However, number of rows and leading dimension are 
             * separated in Chameleon
             */
            int ldauv = BLKLDD(AUV, i);
            int ldad = BLKLDD(AD, i);
            HICMA_TASK_zuncompress(
                    &options,
                    MorseNoTrans, MorseConjTrans,
                    tempi, //number of rows of U 
                    AUV->mb, //number of rows of V 
                    zone,
                    AUV,
                    Ark, //number of columns of U, 
                    // which is equal to that of V
                    i, j, ldauv,
                    zzero,
                    AD, i, j, ldad
                    );
        }
    }
    RUNTIME_sequence_wait( morse, sequence );
    RUNTIME_options_finalize( &options, morse );
    //MORSE_TASK_dataflush_all(); removed in newer chameleon
    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    RUNTIME_options_finalize(&options, morse);




    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}

/**
 * Uncompresses tile low-rank matrix AUV in to AD. 
 * Computes D=U*V^T. Ranks of U and Vs stored in Ark.
 * Number of rows/columns are passed as explicit parameters 
 * to support nonuniform tiles in Tile Low Rank (TLR) format.
 */
int HICMA_zuncompress_custom_size(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark,
        int numrows_matrix, //FIXME use these parameters
        int numcolumns_matrix, //FIXME
        int numrows_block, //FIXME
        int numcolumns_block // FIXME
        )
{

    MORSE_context_t *morse;
    MORSE_sequence_t *sequence = NULL;
    MORSE_request_t request = MORSE_REQUEST_INITIALIZER;
    int status;
    morse = morse_context_self();
    if (morse == NULL) {
        morse_fatal_error("HICMA_zuncompress", "MORSE not initialized");
        return MORSE_ERR_NOT_INITIALIZED;
    }
    morse_sequence_create(morse, &sequence);


    /*MORSE_context_t *morse;*/
    MORSE_option_t options;
    /*morse = morse_context_self();*/
    if (sequence->status != MORSE_SUCCESS)
        return MORSE_ERR_NOT_INITIALIZED;
    RUNTIME_options_init(&options, morse, sequence, &request);
    double zzero = (double) 0.0;
    double zone  = (double) 1.0;
    int i, j;
    for (i = 0; i < AD->mt; i++) {
        /**
         * tempi contains number of rows of tiles in ith row. 
         */
        int tempi = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
        for (j = 0; j < AD->nt; j++) {
            if(uplo == MorseLower && i<=j)
                continue;
            else if(uplo == MorseUpper && i>=j)
                continue;
            int tempj = j == AD->nt-1 ? AD->n-j*AD->nb : AD->nb;
            //printf("%s %d: i:%d j:%d\n", __FILE__,__LINE__,i, j);
            /**
             * leading dimension of AUV. It might be equal to tempi.
             * However, number of rows and leading dimension are 
             * separated in Chameleon
             */
            int ldauv = BLKLDD(AUV, i);
            int ldad = BLKLDD(AD, i);
            HICMA_TASK_zuncompress(
                    &options,
                    MorseNoTrans, MorseConjTrans,
                    tempi, //number of rows of U 
                    tempj, //number of rows of V 
                    zone,
                    AUV,
                    Ark, //number of columns of U, 
                    // which is equal to that of V
                    i, j, ldauv,
                    zzero,
                    AD, i, j, ldad
                    );
        }
    }
    RUNTIME_sequence_wait( morse, sequence );
    RUNTIME_options_finalize( &options, morse );
    //MORSE_TASK_dataflush_all();  removed in newer chameleon
    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    //RUNTIME_options_finalize(&options, morse);

    MORSE_Desc_Flush( AD, sequence );
    MORSE_Desc_Flush( AUV, sequence );
    MORSE_Desc_Flush( Ark, sequence );
    morse_sequence_wait(morse, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

    status = sequence->status;
    morse_sequence_destroy(morse, sequence);
    return status;
}
