/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_async.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU asynchronous routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include "runtime/starpu/hicma_starpu.h"

/**
 *  Create a sequence
 */
int HICMA_RUNTIME_sequence_create( HICMA_context_t  *hicma,
                             HICMA_sequence_t *sequence )
{
    (void)hicma;
    (void)sequence;
    return HICMA_SUCCESS;
}

/**
 *  Destroy a sequence
 */
int HICMA_RUNTIME_sequence_destroy( HICMA_context_t  *hicma,
                              HICMA_sequence_t *sequence )
{
    (void)hicma;
    (void)sequence;
    return HICMA_SUCCESS;
}

/**
 *  Wait for the completion of a sequence
 */
int HICMA_RUNTIME_sequence_wait( HICMA_context_t  *hicma,
                           HICMA_sequence_t *sequence )
{
    (void)hicma;
    (void)sequence;

    if (hicma->progress_enabled) {
        HICMA_RUNTIME_progress(hicma);
    }

    starpu_task_wait_for_all();
#if defined(HICMA_USE_MPI)
    starpu_mpi_barrier(MPI_COMM_WORLD);
#endif
    return HICMA_SUCCESS;
}

/**
 *  Terminate a sequence
 */
void HICMA_RUNTIME_sequence_flush( HICMA_context_t  *hicma,
                             HICMA_sequence_t *sequence,
                             HICMA_request_t  *request,
                             int status )
{
    (void)hicma;
    sequence->request = request;
    sequence->status = status;
    request->status = status;
    return;
}
