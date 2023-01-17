/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file auxiliary.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon auxiliary header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Piotr Luszczek
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#ifndef _HICMA_CHAM_AUXILIARY_H_
#define _HICMA_CHAM_AUXILIARY_H_

#include <hicma_struct.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
void hicma_warning      (const char *func_name, const char* msg_text);
void hicma_error        (const char *func_name, const char* msg_text);
void hicma_fatal_error  (const char *func_name, const char* msg_text);
int  hicma_rank         (HICMA_context_t *hicma);
int  hicma_tune         (HICMA_enum func, int M, int N, int NRHS);

/**
 *  API routines
 */
int  HICMA_Version      (int *ver_major, int *ver_minor, int *ver_micro);
int  HICMA_Element_Size (int type);
int  HICMA_My_Mpi_Rank  (void);

#ifdef __cplusplus
}
#endif

#endif
