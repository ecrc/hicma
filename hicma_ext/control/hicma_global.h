/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file global.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon global variables header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
/**
 *  HICMA internals of interest to HICMA core developers, but not necessarily
 *  of interest to HICMA community contributors.
 */
#ifndef _HICMA_CHAM_GLOBAL_H_
#define _HICMA_CHAM_GLOBAL_H_

#if defined( _WIN32 ) || defined( _WIN64 )
#include "control/hicmawinthread.h"
#else
#include <pthread.h>
#endif

/**
 *  Numerical operations
 */
#define HICMA_FUNC_SGEMM   19
#define HICMA_FUNC_DGEMM   20
#define HICMA_FUNC_CGEMM   21
#define HICMA_FUNC_ZGEMM   22

#endif
