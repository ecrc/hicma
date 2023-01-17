/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file chameleon_config.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon configuration file
 *
 * @version 1.0.0
 * @author Florent Pruvost
 * @date 2017-01-06
 *
 */
#ifndef CHAMELEON_CONFIG_H_HAS_BEEN_INCLUDED
#define CHAMELEON_CONFIG_H_HAS_BEEN_INCLUDED

#define HICMA_CHAM_VERSION_MAJOR 1
#define HICMA_CHAM_VERSION_MINOR 0
#define HICMA_CHAM_VERSION_MICRO 0

/* Scheduling engine */
/* #undef CHAMELEON_SCHED_QUARK */
/* #undef CHAMELEON_SCHED_PARSEC */
#define CHAMELEON_SCHED_STARPU

/* Communication engine */
/* #undef HICMA_USE_MPI */
/* #undef CHAMELEON_USE_MIGRATE */
#if !defined(HICMA_USE_MPI) && defined(CHAMELEON_USE_MIGRATE)
#undef CHAMELEON_USE_MIGRATE
#endif

/* GPU Support */
/* #undef CHAMELEON_USE_CUDA */
/* #undef CHAMELEON_USE_CUBLAS */
/* #undef CHAMELEON_USE_CUBLAS_V2 */

/* Simulation */
/* #undef CHAMELEON_SIMULATION */

/* Tracing support */
/* #undef CHAMELEON_ENABLE_TRACING */

/* getopt */
#define CHAMELEON_HAVE_GETOPT_H
#define CHAMELEON_HAVE_GETOPT_LONG

#ifdef BEGIN_C_DECLS
#undef BEGIN_C_DECLS
#endif

#ifdef END_C_DECLS
#undef END_C_DECLS
#endif

#if defined(c_plusplus) || defined(__cplusplus)
# define BEGIN_C_DECLS extern "C" {
# define END_C_DECLS }
#else
# define BEGIN_C_DECLS          /* empty */
# define END_C_DECLS            /* empty */
#endif

#endif  /* CHAMELEON_CONFIG_H_HAS_BEEN_INCLUDED */
