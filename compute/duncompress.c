/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 *
 * @file duncompress.c
 *
 * This file contains the function for uncompressing a tile low-rank matrix.
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Kadir Akbudak
 * @date 2018-11-08
 *
 **/

#include <control/common.h>
#include <include/hicma_d.h>
#include <include/hicma_runtime_d.h>

/**
 * Uncompresses tile low-rank matrix AUV in to AD. 
 * Computes D=U*V^T. Ranks of U and Vs stored in Ark
 */
int HICMA_duncompress(HICMA_enum uplo,
        HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark)
{

    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;
    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_duncompress", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);


    /*HICMA_context_t *hicma;*/
    HICMA_option_t options;
    /*hicma = hicma_context_self();*/
    if (sequence->status != HICMA_SUCCESS)
        return HICMA_ERR_NOT_INITIALIZED;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, &request);
    double dzero = (double) 0.0;
    double done  = (double) 1.0;
    int i, j;
    for (i = 0; i < AD->mt; i++) {
        /**
         * tempi contains number of rows of tiles in ith row. 
         */
        int tempi = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
        for (j = 0; j < AD->mt; j++) {
            if(uplo == HicmaLower && i<=j)
                continue;
            else if(uplo == HicmaUpper && i>=j)
                continue;
            //printf("%s %d: i:%d j:%d\n", __FILE__,__LINE__,i, j);
            /**
             * leading dimension of AUV. It might be equal to tempi.
             * However, number of rows and leading dimension are 
             * separated in Chameleon
             */
            int ldauv = BLKLDD(AUV, i);
            int ldad = BLKLDD(AD, i);
            HICMA_TASK_duncompress(
                    &options,
                    HicmaNoTrans, HicmaConjTrans,
                    tempi, //number of rows of U 
                    AUV->mb, //number of rows of V 
                    done,
                    AUV,
                    Ark, //number of columns of U, 
                    // which is equal to that of V
                    i, j, ldauv,
                    dzero,
                    AD, i, j, ldad
                    );
        }
    }
    HICMA_RUNTIME_sequence_wait( hicma, sequence );
    HICMA_RUNTIME_options_finalize( &options, hicma );
    //HICMA_TASK_dataflush_all(); removed in newer chameleon
    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    HICMA_RUNTIME_options_finalize(&options, hicma);




    HICMA_Desc_Flush( AD, sequence );
    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}

/**
 * Uncompresses tile low-rank matrix AUV in to AD. 
 * Computes D=U*V^T. Ranks of U and Vs stored in Ark.
 * Number of rows/columns are passed as explicit parameters 
 * to support nonuniform tiles in Tile Low Rank (TLR) format.
 */
int HICMA_duncompress_custom_size(HICMA_enum uplo,
        HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark,
        int numrows_matrix, //FIXME use these parameters
        int numcolumns_matrix, //FIXME
        int numrows_block, //FIXME
        int numcolumns_block // FIXME
        )
{

    HICMA_context_t *hicma;
    HICMA_sequence_t *sequence = NULL;
    HICMA_request_t request = HICMA_REQUEST_INITIALIZER;
    int status;
    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_duncompress", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    hicma_sequence_create(hicma, &sequence);


    /*HICMA_context_t *hicma;*/
    HICMA_option_t options;
    /*hicma = hicma_context_self();*/
    if (sequence->status != HICMA_SUCCESS)
        return HICMA_ERR_NOT_INITIALIZED;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, &request);
    double dzero = (double) 0.0;
    double done  = (double) 1.0;
    int i, j;
    for (i = 0; i < AD->mt; i++) {
        /**
         * tempi contains number of rows of tiles in ith row. 
         */
        int tempi = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
        for (j = 0; j < AD->nt; j++) {
            if(uplo == HicmaLower && i<=j)
                continue;
            else if(uplo == HicmaUpper && i>=j)
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
            HICMA_TASK_duncompress(
                    &options,
                    HicmaNoTrans, HicmaConjTrans,
                    tempi, //number of rows of U 
                    tempj, //number of rows of V 
                    done,
                    AUV,
                    Ark, //number of columns of U, 
                    // which is equal to that of V
                    i, j, ldauv,
                    dzero,
                    AD, i, j, ldad
                    );
        }
    }
    HICMA_RUNTIME_sequence_wait( hicma, sequence );
    HICMA_RUNTIME_options_finalize( &options, hicma );
    //HICMA_TASK_dataflush_all();  removed in newer chameleon
    //RUNTIME_desc_getoncpu( &AD ); accuracy checking works without this line on shared memory and with 4 mpi ranks on shared memory
    //HICMA_RUNTIME_options_finalize(&options, hicma);

    HICMA_Desc_Flush( AD, sequence );
    HICMA_Desc_Flush( AUV, sequence );
    HICMA_Desc_Flush( Ark, sequence );
    hicma_sequence_wait(hicma, sequence);
    /*RUNTIME_desc_getoncpu(AD);*/
    /*RUNTIME_desc_getoncpu(AUV);*/
    /*RUNTIME_desc_getoncpu(Ark);*/

    status = sequence->status;
    hicma_sequence_destroy(hicma, sequence);
    return status;
}
