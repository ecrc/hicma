/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file hicma_types.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon basic datatypes header
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _HICMA_CHAM_TYPES_H_
#define _HICMA_CHAM_TYPES_H_

#include "hicma_config.h"

/**
 *  System requirements
 */
#include <stddef.h>
#if defined( _WIN32 )
  /* This must be included before INPUT is defined below, otherwise we
     have a name clash/problem  */
  #include <windows.h>
  #include <limits.h>
#else /* _WIN32 */
  #include <inttypes.h>
#endif /* _WIN32 */


/**
 *  HICMA types
 */
typedef int  HICMA_enum;
typedef int  HICMA_bool;
typedef long HICMA_index;
typedef long HICMA_size;


/**
 * HICMA Complex numbers
 */
#define HICMA_HAS_COMPLEX_H 1

#if defined(_WIN32)
# include <float.h>
# if defined(__INTEL_COMPILER)
    /* Fix name conflict within the cabs prototype (_Complex) that    */
    /* conflicts with a C99 keyword.                                  */
    #define _Complex __ConflictingComplex
    #include <math.h>
    #undef _Complex
    #undef complex
# elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    #undef  HICMA_COMPLEX_CPP
    #define HICMA_COMPLEX_CPP
# else
    #error "Supported compilers on WIN32 are MSVC and Intel Compiler."
# endif /* __INTEL_COMPILER */

# define isnan _isnan
# define isinf !_finite
#endif /* _WIN32 */

/* Sun doesn't ship the complex.h header. Sun Studio doesn't have it and older GCC compilers don't have it either. */
#if defined(__SUNPRO_C) || defined(__SUNPRO_CC) || defined(sun) || defined(__sun)
#undef HICMA_HAS_COMPLEX_H
#endif /* __SUNPRO_C */

#ifndef __cplusplus
    #undef HICMA_COMPLEX_CPP
#endif

#if defined(HICMA_COMPLEX_CPP)
    #ifndef LAPACK_COMPLEX_CPP
    # define LAPACK_COMPLEX_CPP
    # warning "HiCMA_COMPLEX_CPP was defined, but not LAPACK_COMPLEX_CPP. Maybe you want to set both."
    #endif
    #include <complex> // needed for std::complex declaration
    #define HICMA_Complex32_t std::complex<float>
    #define HICMA_Complex64_t std::complex<double>
#else /* HICMA_COMPLEX_CPP */
      /* not using cplusplus complex type: */

    #if defined(__STDC_NO_COMPLEX__)
    # error "Compiler support for complex number is required."
    #endif

    #define HICMA_Complex32_t float  _Complex
    #define HICMA_Complex64_t double _Complex

    #if HICMA_HAS_COMPLEX_H
    # include <complex.h>
    #endif
#endif /* HICMA_COMPLEX_CPP */

/**
 *  HICMA Deprecated attribute
 */
#if (__GNUC__ > 3 || (__GNUC__ == 3 && __GNUC_MINOR__ >= 1))
#define HICMA_DEPRECATED  __attribute__((__deprecated__))
#else
#define HICMA_DEPRECATED
#endif /* __GNUC__ */

BEGIN_C_DECLS

/**
 *  Global utilities
 */
static inline int hicma_max( int a, int b ) {
    if ( a > b ) return a; else return b;
}

static inline int hicma_min( int a, int b ) {
    if ( a < b ) return a; else return b;
}

END_C_DECLS

#endif /* __CHAMELEON_H__ */
