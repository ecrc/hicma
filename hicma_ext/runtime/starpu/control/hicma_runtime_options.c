/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_options.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU options routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include "runtime/starpu/hicma_starpu.h"

void HICMA_RUNTIME_options_init( HICMA_option_t *option, HICMA_context_t *hicma,
                           HICMA_sequence_t *sequence, HICMA_request_t *request )
{
    option->sequence   = sequence;
    option->request    = request;
    option->profiling  = HICMA_PROFILING == HICMA_TRUE;
    option->parallel   = HICMA_PARALLEL == HICMA_TRUE;
    option->priority   = HICMA_PRIORITY_MIN;
    option->nb         = HICMA_NB;
    option->ws_wsize   = 0;
    option->ws_hsize   = 0;
    option->ws_worker  = NULL;
    option->ws_host    = NULL;
    return;
}

void HICMA_RUNTIME_options_finalize( HICMA_option_t *option, HICMA_context_t *hicma )
{
    (void)option;
    (void)hicma;
    return;
}

int HICMA_RUNTIME_options_ws_alloc( HICMA_option_t *options, size_t worker_size, size_t host_size )
{
    int ret = 0;
    if ( worker_size > 0 ) {
        options->ws_wsize = worker_size;
        starpu_vector_data_register((starpu_data_handle_t*)(&(options->ws_worker)),
                                    -1, (uintptr_t)NULL,
                                    worker_size, sizeof(char));
    }
    if ( host_size > 0 ) {
        options->ws_hsize = host_size;
        ret = HICMA_RUNTIME_starpu_ws_alloc((HICMA_starpu_ws_t**)&(options->ws_host),
                                      host_size, HICMA_CUDA, HICMA_HOST_MEM);
    }
    return ret;
}

int HICMA_RUNTIME_options_ws_free( HICMA_option_t *options )
{
    int ret = 0;
    if ( options->ws_worker != NULL ) {
        starpu_data_unregister_submit((starpu_data_handle_t)(options->ws_worker));
        options->ws_worker = NULL;
    }
    if ( options->ws_host != NULL ) {
        starpu_task_wait_for_all();
        ret = HICMA_RUNTIME_starpu_ws_free( (HICMA_starpu_ws_t*)(options->ws_host) );
        options->ws_host = NULL;
    }
    return ret;
}
