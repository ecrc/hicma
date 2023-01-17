/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Sequences
 * @brief Group routines exposed to users to handle asynchronous tasks execution
 *
 */
#include <stdlib.h>
#include "common.h"
#include <hicma_runtime.h>

/**
 *  Register an exception.
 */
int hicma_request_fail(HICMA_sequence_t *sequence, HICMA_request_t *request, int status)
{
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return status;
}

/**
 *  Create a sequence
 */
int hicma_sequence_create(HICMA_context_t *hicma, HICMA_sequence_t **sequence)
{
    if ((*sequence = malloc(sizeof(HICMA_sequence_t))) == NULL) {
        hicma_error("HiCMA_Sequence_Create", "malloc() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }

    HICMA_RUNTIME_sequence_create( hicma, *sequence );

    (*sequence)->status = HICMA_SUCCESS;
    return HICMA_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int hicma_sequence_destroy(HICMA_context_t *hicma, HICMA_sequence_t *sequence)
{
    HICMA_RUNTIME_sequence_destroy( hicma, sequence );
    free(sequence);
    return HICMA_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int hicma_sequence_wait(HICMA_context_t *hicma, HICMA_sequence_t *sequence)
{
    HICMA_RUNTIME_sequence_wait( hicma, sequence );
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Sequences
 *
 *  HICMA_Sequence_Create - Create a squence.
 *
 ******************************************************************************
 *
 * @param[out] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Sequence_Create(HICMA_sequence_t **sequence)
{
    HICMA_context_t *hicma;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Create", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    status = hicma_sequence_create(hicma, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  HICMA_Sequence_Destroy - Destroy a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Sequence_Destroy(HICMA_sequence_t *sequence)
{
    HICMA_context_t *hicma;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Destroy", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Destroy", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    status = hicma_sequence_destroy(hicma, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  HICMA_Sequence_Wait - Wait for the completion of a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Sequence_Wait(HICMA_sequence_t *sequence)
{
    HICMA_context_t *hicma;
    int status;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Wait", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Wait", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }
    status = hicma_sequence_wait(hicma, sequence);
    return status;
}

/**
 *
 * @ingroup Sequences
 *
 *  HICMA_Sequence_Flush - Terminate a sequence.
 *
 ******************************************************************************
 *
 * @param[in] sequence
 *          Identifies a set of routines sharing common exception handling.
 *
 * @param[in] request
 *          The flush request.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Sequence_Flush(HICMA_sequence_t *sequence, HICMA_request_t *request)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Flush", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (sequence == NULL) {
        hicma_fatal_error("HiCMA_Sequence_Flush", "NULL sequence");
        return HICMA_ERR_UNALLOCATED;
    }

    HICMA_RUNTIME_sequence_flush( hicma->schedopt, sequence, request, HICMA_ERR_SEQUENCE_FLUSHED);

    return HICMA_SUCCESS;
}
