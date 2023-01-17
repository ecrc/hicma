/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file coreblas.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon CPU kernels main header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @date 2010-11-15
 *
 */
#ifndef _CORE_BLAS_H_
#define _CORE_BLAS_H_

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>

/**
 *  CBLAS requires for scalar arguments to be passed
 *        by address rather than by value
 */
#ifndef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
#endif
#include "cblas.h"

/**
 * HICMA types and constants
 */
#include <hicma_config.h>
#include <hicma_types.h>
#include <hicma_struct.h>
#include "hicma_constants.h"


/**
 * Coreblas Error
 */
#define coreblas_error(k, str) do {                                     \
        fprintf(stderr, "%s: Parameter %d / %s\n", __func__, k, str) ;  \
        assert(0);                                                      \
    } while(0)

/**
 * CBlas enum
 */
#define CBLAS_TRANSPOSE enum CBLAS_TRANSPOSE
#define CBLAS_UPLO      enum CBLAS_UPLO
#define CBLAS_DIAG      enum CBLAS_DIAG
#define CBLAS_SIDE      enum CBLAS_SIDE

/**
 *  LAPACK Constants
 */
BEGIN_C_DECLS

extern char *hicma_lapack_constants[];
#define hicma_lapack_const(hicma_const) hicma_lapack_constants[hicma_const][0]

void HICMA_set_coreblas_gemm3m_enabled( int v );
int  HICMA_get_coreblas_gemm3m_enabled( void );

END_C_DECLS

#endif
