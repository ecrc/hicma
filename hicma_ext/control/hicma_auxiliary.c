/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file auxiliary.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary routines
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 ***
 *
 * @defgroup Auxiliary
 * @brief Group auxiliary routines exposed to users
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include "common.h"
#include "hicma_auxiliary.h"

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action without severe consequences.
 *  Problems occuring while HICMA is being used correctly.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void hicma_warning(const char *func_name, const char *msg_text)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL)
        hicma_fatal_error("hicma_warning", "HiCMA not initialized");
    if (hicma->warnings_enabled)
        fprintf(stderr, "HiCMA WARNING: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Indicates a recoverable problem.
 *  User's erroneous action with potentially severe consequences.
 *  Problems occuring due to incorrect use of HICMA.
 *  Context aware.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void hicma_error(const char *func_name, const char *msg_text)
{
    fprintf(stderr, "HiCMA ERROR: %s(): %s\n", func_name, msg_text);
}

/**
 *
 *  Unexpected behavior within the library.
 *  Unrecoverable user errors.
 *  Context oblivious.
 *
 * @param[in] func_name
 *          Function location where warning occurred
 *
 * @param[in] msg_text
 *          Warning message to display.
 *
 */
void hicma_fatal_error(const char *func_name, const char *msg_text)
{
    fprintf(stderr, "HiCMA FATAL ERROR: %s(): %s\n", func_name, msg_text);
    exit(0);
}

/**
 *  Returns core id
 */
int hicma_rank(HICMA_context_t *hicma)
{
    return HICMA_RUNTIME_thread_rank( hicma );
}

/**
 *  Tune block size nb and internal block size ib
 */
int hicma_tune(HICMA_enum func, int M, int N, int NRHS)
{
    HICMA_context_t *hicma;
    hicma = hicma_context_self();
    if ( hicma && hicma->autotuning_enabled == HICMA_TRUE ) {
        hicma_warning( "hicma_tune", "Autotunning not available for now" );
    }
    (void)func;
    (void)M;
    (void)N;
    (void)NRHS;
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  HICMA_Version - Reports HICMA version number.
 *
 ******************************************************************************
 *
 * @param[out] ver_major
 *          HICMA major version number.
 *
 * @param[out] ver_minor
 *          HICMA minor version number.
 *
 * @param[out] ver_micro
 *          HICMA micro version number.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Version(int *ver_major, int *ver_minor, int *ver_micro)
{
    if (! ver_major && ! ver_minor && ! ver_micro)
        return  HICMA_ERR_ILLEGAL_VALUE;

    if (ver_major)
        *ver_major = HICMA_CHAM_VERSION_MAJOR;

    if (ver_minor)
        *ver_minor = HICMA_CHAM_VERSION_MINOR;

    if (ver_micro)
        *ver_micro = HICMA_CHAM_VERSION_MICRO;

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Auxiliary
 *
 *  HICMA_Element_Size - Reports the size in bytes of a HICMA precision type
 *  (e.g. HicmaInteger, HicmaRealFloat, etc).
 *
 ******************************************************************************
 *
 * @param[in] type
 *          HICMA element type, can be one of the following:
 *          - HicmaByte
 *          - HicmaInteger
 *          - HicmaRealFloat
 *          - HicmaRealDouble
 *          - HicmaComplexFloat
 *          - HicmaComplexDouble
 *
 ******************************************************************************
 *
 * @return
 *          \retval Element size in bytes
 *
 */
int HICMA_Element_Size(int type)
{
    switch(type) {
        case HicmaByte:          return          1;
        case HicmaInteger:       return   sizeof(int);
        case HicmaRealFloat:     return   sizeof(float);
        case HicmaRealDouble:    return   sizeof(double);
        case HicmaComplexFloat:  return 2*sizeof(float);
        case HicmaComplexDouble: return 2*sizeof(double);
        default: hicma_fatal_error("HiCMA_Element_Size", "undefined type");
                 return HICMA_ERR_ILLEGAL_VALUE;

    }
}

/**
 *
 * @ingroup Auxiliary
 *
 *  HICMA_My_Mpi_Rank - Return the MPI rank of the calling process.
 *
 ******************************************************************************
 *
 ******************************************************************************
 *
 * @return
 *          \retval MPI rank
 *
 */
int HICMA_My_Mpi_Rank(void)
{
#if defined(HICMA_USE_MPI)
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Finalize()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    return HICMA_MPI_RANK;
#else
    return HICMA_SUCCESS;
#endif
}

/**
 *  Display a progress percentage in stderr
 */
void HICMA_update_progress(int currentValue, int maximumValue) {
    div_t res ;
    static int progress = -1; /* varie de 0 a 100 au cours du calcul concerne */

    if (maximumValue == 0) {
        res.quot = 100;
    }
    else {
        if (currentValue < (INT_MAX / 100) ) {
            res = div(currentValue*100, maximumValue);
        }
        else {
            /* Calcule le quotient de la division */
            res.quot = (int)( (long long)( currentValue * 100 ) / maximumValue );
        }
    }

    // Print the percentage
    if (res.quot > progress) {
        fprintf(stderr, "%3d%%\b\b\b\b", res.quot);
    }
    progress = res.quot;

    if (currentValue >= maximumValue) {
        progress = -1;
    }
}

// A function to display the progress indicator.
// By default it is HICMA_update_progress()
// The user can change it with HICMA_Set_HICMA_update_progress_callback()
void (*HICMA_update_progress_callback)(int, int) = HICMA_update_progress;

int HICMA_Set_HICMA_update_progress_callback(void (*p)(int, int)) {
  HICMA_update_progress_callback = p;
  return HICMA_SUCCESS;
}

