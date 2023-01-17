/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file context.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context management routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Options
 * @brief Group routines exposed to users to handle options
 *
 */

#include <stdlib.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include "control/morsewinthread.h"
#else
#include <pthread.h>
#endif

#include "common.h"
#include "hicma_auxiliary.h"
#include "hicma_context.h"
#include <hicma_runtime.h>

#if !defined(CHAMELEON_SIMULATION)
#include "coreblas/hicma_coreblas.h"
#endif

/**
 *  Global data
 */
/* master threads context lookup table */
static HICMA_context_t *hicma_ctxt = NULL;

/**
 *  Create new context
 */
HICMA_context_t *hicma_context_create()
{
    HICMA_context_t *hicma;

    if ( hicma_ctxt != NULL ) {
        hicma_error("hicma_context_create", "a context is already existing\n");
        return NULL;
    }

    hicma = (HICMA_context_t*)malloc(sizeof(HICMA_context_t));
    if (hicma == NULL) {
        hicma_error("hicma_context_create", "malloc() failed");
        return NULL;
    }

    /* These initializations are just in case the user
       disables autotuning and does not set nb and ib */
    hicma->nb                 = 128;
    hicma->ib                 = 32;
    hicma->rhblock            = 4;

    hicma->nworkers           = 1;
    hicma->ncudas             = 0;
    hicma->nthreads_per_worker= 1;

    hicma->warnings_enabled     = HICMA_TRUE;
    hicma->autotuning_enabled   = HICMA_TRUE;
    hicma->parallel_enabled     = HICMA_FALSE;
    hicma->profiling_enabled    = HICMA_FALSE;
    hicma->progress_enabled     = HICMA_FALSE;

    hicma->householder        = HICMA_FLAT_HOUSEHOLDER;
    hicma->translation        = HICMA_OUTOFPLACE;


    /* Initialize scheduler */
    HICMA_RUNTIME_context_create(hicma);

    hicma_ctxt = hicma;
    return hicma;
}


/**
 *  Return context for a thread
 */
HICMA_context_t *hicma_context_self()
{
    return hicma_ctxt;
}

/**
 *  Clean the context
 */
int hicma_context_destroy(){

    HICMA_RUNTIME_context_destroy(hicma_ctxt);
    free(hicma_ctxt);
    hicma_ctxt = NULL;

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  HICMA_Enable - Enable HICMA feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be enabled:
 *          @arg HICMA_WARNINGS   printing of warning messages,
 *          @arg HICMA_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg HICMA_PROFILING_MODE  activate profiling of kernels
 *          @arg HICMA_PROGRESS  activate progress indicator
 *          @arg HICMA_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Enable(HICMA_enum option)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Enable", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }

    switch (option)
    {
        case HICMA_WARNINGS:
            hicma->warnings_enabled = HICMA_TRUE;
            break;
        case HICMA_AUTOTUNING:
            hicma->autotuning_enabled = HICMA_TRUE;
            break;
        case HICMA_PROFILING_MODE:
            hicma->profiling_enabled = HICMA_TRUE;
            break;
        case HICMA_PROGRESS:
            hicma->progress_enabled = HICMA_TRUE;
            break;
        case HICMA_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            HICMA_set_coreblas_gemm3m_enabled(1);
#else
            hicma_error("HiCMA_Enable", "cannot enable GEMM3M (not available in cblas)");
#endif
            break;
        /* case HICMA_PARALLEL: */
        /*     hicma->parallel_enabled = HICMA_TRUE; */
        /*     break; */
        default:
            hicma_error("HiCMA_Enable", "illegal parameter value");
            return HICMA_ERR_ILLEGAL_VALUE;
        case HICMA_BOUND:
            break;
    }

    /* Enable at the lower level if required */
    HICMA_RUNTIME_enable( option );

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  HICMA_Disable - Disable HICMA feature.
 *
 *******************************************************************************
 *
 * @param[in] option
 *          Feature to be disabled:
 *          @arg HICMA_WARNINGS   printing of warning messages,
 *          @arg HICMA_AUTOTUNING autotuning for tile size and inner block size.
 *          @arg HICMA_PROFILING_MODE  deactivate profiling of kernels
 *          @arg HICMA_PROGRESS  deactivate progress indicator
 *          @arg HICMA_GEMM3M  Use z/cgemm3m for complexe matrix-matrix products
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Disable(HICMA_enum option)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Disable", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    switch ( option )
    {
        case HICMA_WARNINGS:
            hicma->warnings_enabled = HICMA_FALSE;
            break;
        case HICMA_AUTOTUNING:
            hicma->autotuning_enabled = HICMA_FALSE;
            break;
        case HICMA_PROFILING_MODE:
            hicma->profiling_enabled = HICMA_FALSE;
            break;
        case HICMA_PROGRESS:
            hicma->progress_enabled = HICMA_FALSE;
            break;
        case HICMA_GEMM3M:
#if defined(CBLAS_HAS_ZGEMM3M) && !defined(CHAMELEON_SIMULATION)
            HICMA_set_coreblas_gemm3m_enabled(0);
#endif
            break;
        case HICMA_PARALLEL_MODE:
            hicma->parallel_enabled = HICMA_FALSE;
            break;
        default:
            hicma_error("HiCMA_Disable", "illegal parameter value");
            return HICMA_ERR_ILLEGAL_VALUE;
    }

    /* Disable at the lower level if required */
    HICMA_RUNTIME_disable( option );

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  HICMA_Set - Set HICMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg HICMA_TILE_SIZE:        size matrix tile,
 *          @arg HICMA_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[in] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Set(HICMA_enum param, int value)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Set", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case HICMA_TILE_SIZE:
            if (value <= 0) {
                hicma_error("HiCMA_Set", "negative tile size");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            hicma->nb = value;
            if ( hicma->autotuning_enabled ) {
                hicma->autotuning_enabled = HICMA_FALSE;
                hicma_warning("HiCMA_Set", "autotuning has been automatically disable\n");
            }
            /* Limit ib to nb */
            hicma->ib = hicma_min( hicma->nb, hicma->ib );
            break;
        case HICMA_INNER_BLOCK_SIZE:
            if (value <= 0) {
                hicma_error("HiCMA_Set", "negative inner block size");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            if (value > hicma->nb) {
                hicma_error("HiCMA_Set", "inner block larger than tile");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            /* if (hicma->nb % value != 0) { */
            /*     hicma_error("HiCMA_Set", "inner block does not divide tile"); */
            /*     return HICMA_ERR_ILLEGAL_VALUE; */
            /* } */
            hicma->ib = value;

            if ( hicma->autotuning_enabled ) {
                hicma->autotuning_enabled = HICMA_FALSE;
                hicma_warning("HiCMA_Set", "autotuning has been automatically disable\n");
            }
            break;
        case HICMA_HOUSEHOLDER_MODE:
            if (value != HICMA_FLAT_HOUSEHOLDER && value != HICMA_TREE_HOUSEHOLDER) {
                hicma_error("HiCMA_Set", "illegal value of HICMA_HOUSEHOLDER_MODE");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            hicma->householder = value;
            break;
        case HICMA_HOUSEHOLDER_SIZE:
            if (value <= 0) {
                hicma_error("HiCMA_Set", "negative householder size");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            hicma->rhblock = value;
            break;
        case HICMA_TRANSLATION_MODE:
            if (value != HICMA_INPLACE && value != HICMA_OUTOFPLACE) {
                hicma_error("HiCMA_Set", "illegal value of HICMA_TRANSLATION_MODE");
                return HICMA_ERR_ILLEGAL_VALUE;
            }
            hicma->translation = value;
            break;
        default:
            hicma_error("HiCMA_Set", "unknown parameter");
            return HICMA_ERR_ILLEGAL_VALUE;
    }

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Options
 *
 *  HICMA_Get - Get value of HICMA parameter.
 *
 *******************************************************************************
 *
 * @param[in] param
 *          Feature to be enabled:
 *          @arg HICMA_TILE_SIZE:        size matrix tile,
 *          @arg HICMA_INNER_BLOCK_SIZE: size of tile inner block,
 *
 * @param[out] value
 *          Value of the parameter.
 *
 *******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Get(HICMA_enum param, int *value)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Get", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    switch (param) {
        case HICMA_TILE_SIZE:
            *value = hicma->nb;
            return HICMA_SUCCESS;
        case HICMA_INNER_BLOCK_SIZE:
            *value = hicma->ib;
            return HICMA_SUCCESS;
        case HICMA_HOUSEHOLDER_MODE:
            *value = hicma->householder;
            return HICMA_SUCCESS;
        case HICMA_HOUSEHOLDER_SIZE:
            *value = hicma->rhblock;
            return HICMA_SUCCESS;
        case HICMA_TRANSLATION_MODE:
            *value = hicma->translation;
            return HICMA_SUCCESS;
        default:
            hicma_error("HiCMA_Get", "unknown parameter");
            return HICMA_ERR_ILLEGAL_VALUE;
    }

    return HICMA_SUCCESS;
}
