/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU context routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include "runtime/starpu/hicma_starpu.h"

#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 3))
/* Defined by StarPU as external function */
#else
#if ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION >= 2))
int _starpu_is_initialized(void);
#define starpu_is_initialized() _starpu_is_initialized()
#else
#define starpu_is_initialized() 0
#endif
#endif

/**
 *  Create new context
 */
void HICMA_RUNTIME_context_create( HICMA_context_t *hicma )
{
    starpu_conf_t *conf;

    hicma->scheduler = HICMA_RUNTIME_SCHED_STARPU;

    if (! starpu_is_initialized() ) {
        hicma->schedopt = (void*) malloc (sizeof(starpu_conf_t));
        conf = hicma->schedopt;

        starpu_conf_init( conf );
    }
    else {
        hicma->schedopt = NULL;
    }

    return;
}

/**
 *  Clean the context
 */
void HICMA_RUNTIME_context_destroy( HICMA_context_t *hicma )
{
    /* StarPU was already initialized by an external library */
    if (hicma->schedopt) {
        free(hicma->schedopt);
    }
    return;
}

/**
 *
 */
void HICMA_RUNTIME_enable( HICMA_enum lever )
{
    switch (lever)
    {
        case HICMA_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_ENABLE);
            break;
        case HICMA_BOUND:
            starpu_bound_start(0, 0);
            break;
        default:
            return;
    }
    return;
}

/**
 *
 */
void HICMA_RUNTIME_disable( HICMA_enum lever )
{
    switch (lever)
    {
        case HICMA_PROFILING_MODE:
            starpu_profiling_status_set(STARPU_PROFILING_DISABLE);
            break;
        case HICMA_BOUND:
            starpu_bound_stop();
            break;
        default:
            return;
    }
    return;
}
