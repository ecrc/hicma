/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file async.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon asynchronous management header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#ifndef _HICMA_CHAM_ASYNC_H_
#define _HICMA_CHAM_ASYNC_H_

#include <hicma_struct.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
int hicma_request_fail     (HICMA_sequence_t *sequence, HICMA_request_t *request, int error);
int hicma_sequence_create  (HICMA_context_t *HICMA, HICMA_sequence_t **sequence);
int hicma_sequence_destroy (HICMA_context_t *HICMA, HICMA_sequence_t *sequence);
int hicma_sequence_wait    (HICMA_context_t *HICMA, HICMA_sequence_t *sequence);

#ifdef __cplusplus
}
#endif

#endif
