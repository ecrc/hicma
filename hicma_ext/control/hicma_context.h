/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file context.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon context header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#ifndef _HICMA_CHAM_CONTEXT_H_
#define _HICMA_CHAM_CONTEXT_H_

#include <hicma_struct.h>

/**
 *  Routines to handle threads context
 */
#ifdef __cplusplus
extern "C" {
#endif

HICMA_context_t* hicma_context_create  ();
HICMA_context_t* hicma_context_self    ();
int              hicma_context_destroy ();

#ifdef __cplusplus
}
#endif

#endif
