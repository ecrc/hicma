/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file chameleon_starpu.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU runtime header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2011-06-01
 *
 */
#ifndef _HICMA_CHAM_STARPU_H_
#define _HICMA_CHAM_STARPU_H_

#include <hicma_config.h>

/* StarPU options */
/* #undef HAVE_STARPU_FXT_PROFILING */
/* #undef HAVE_STARPU_IDLE_PREFETCH */
/* #undef HAVE_STARPU_ITERATION_PUSH */
/* #undef HAVE_STARPU_DATA_WONT_USE */
/* #undef HAVE_STARPU_DATA_SET_COORDINATES */
/* #undef HAVE_STARPU_MALLOC_ON_NODE_SET_DEFAULT_FLAGS */
/* #undef HAVE_STARPU_MPI_DATA_MIGRATE */
/* #undef HAVE_STARPU_MPI_DATA_REGISTER */
/* #undef HAVE_STARPU_MPI_COMM_RANK */
/* #undef HAVE_STARPU_MPI_CACHED_RECEIVE */
/* #undef HAVE_STARPU_MPI_COMM_GET_ATTR */

#if defined(HICMA_USE_MPI)
#include <starpu_mpi.h>
#else
#include <starpu.h>
#endif

#include <starpu_profiling.h>

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <starpu_scheduler.h>
#include <starpu_cuda.h>

#include <cublas.h>
#include <starpu_cublas.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas_v2.h>
#include <starpu_cublas_v2.h>
#endif
#endif

#if defined(CHAMELEON_SIMULATION)
# if !defined(STARPU_SIMGRID)
#  error "Starpu was not built with simgrid support (--enable-simgrid). Can not run Chameleon with simulation support."
# endif
#else
# if defined(STARPU_SIMGRID)
#  warning "Starpu was built with simgrid support. Better build Chameleon with simulation support (-DCHAMELEON_SIMULATION=YES)."
# endif
#endif

#include <control/common.h>
#include "runtime/starpu/hicma_runtime_codelets.h"
#include "runtime/starpu/hicma_runtime_profiling.h"
#include "runtime/starpu/hicma_runtime_codelet_profile.h"
#include "runtime/starpu/hicma_runtime_workspace.h"

typedef struct starpu_conf starpu_conf_t;

/**/

/*
 * MPI Redefinitions
 */
#if defined(HICMA_USE_MPI)
#undef STARPU_REDUX
//#define starpu_insert_task(...) starpu_mpi_insert_task(MPI_COMM_WORLD, __VA_ARGS__)
#define starpu_insert_task starpu_mpi_insert_task
#define starpu_mpi_codelet(_codelet_) MPI_COMM_WORLD, _codelet_

#else

#define starpu_mpi_codelet(_codelet_) _codelet_

#endif

/*
 * cuBlasAPI v2 - StarPU enable the support for cublas handle
 */
#if defined(CHAMELEON_USE_CUDA) && defined(CHAMELEON_USE_CUBLAS_V2)
#define RUNTIME_getStream(_stream_)                             \
    cublasHandle_t _stream_ = starpu_cublas_get_local_handle();
#else
#define RUNTIME_getStream(_stream_)                             \
    cudaStream_t _stream_ = starpu_cuda_get_local_stream();     \
    cublasSetKernelStream( stream );

#endif

/*
 * Enable codelets names
 */
#if (STARPU_MAJOR_VERSION > 1) || ((STARPU_MAJOR_VERSION == 1) && (STARPU_MINOR_VERSION > 1))
#define CHAMELEON_CODELETS_HAVE_NAME
#endif

/**
 * Access to block pointer and leading dimension
 */
#define RTBLKADDR( desc, type, m, n ) ( (starpu_data_handle_t)HICMA_RUNTIME_data_getaddr( desc, m, n ) )

void RUNTIME_set_reduction_methods(starpu_data_handle_t handle, HICMA_enum dtyp);

#if defined(HICMA_USE_MPI) && defined(HAVE_STARPU_MPI_CACHED_RECEIVE)
static inline int
chameleon_starpu_data_iscached(const HICMA_desc_t *A, int m, int n)
{
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t)A->lmt) * nn + mm;

    if (!(*ptrtile))
        return 0;

    return starpu_mpi_cached_receive(*ptrtile);
}

#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {                 \
        if (chameleon_starpu_data_iscached(A, Am, An)) __hicma_need_submit = 1; } while(0)

#else

#warning "WAR dependencies need starpu_mpi_cached_receive support from StarPU 1.2.1 or greater"
#define RUNTIME_ACCESS_WRITE_CACHED(A, Am, An) do {} while (0)

#endif

#ifdef CHAMELEON_ENABLE_PRUNING_STATS

#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION \
    int __hicma_exec = 0; \
    int __hicma_changed = 0;

#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An) \
    if (hicma_desc_islocal(A, Am, An)) \
        __hicma_exec = 1;

#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION \
    RUNTIME_total_tasks++; \
    if (__hicma_exec) \
        RUNTIME_exec_tasks++; \
    else if (__hicma_need_submit) \
        RUNTIME_comm_tasks++; \
    else if (__hicma_changed) \
        RUNTIME_changed_tasks++;

#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank) \
    int __hicma_myrank; \
    HICMA_RUNTIME_comm_rank(&__hicma_myrank); \
    __hicma_exec = (rank) == __hicma_myrank; \
    __hicma_changed = 1; \

#else
#define RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An)
#define RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION
#define RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)
#endif

#define RUNTIME_BEGIN_ACCESS_DECLARATION        \
    RUNTIME_PRUNING_STATS_BEGIN_ACCESS_DECLARATION

#define RUNTIME_ACCESS_R(A, Am, An)

#define RUNTIME_ACCESS_W(A, Am, An)             \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_ACCESS_RW(A, Am, An)            \
    RUNTIME_PRUNING_STATS_ACCESS_W(A, Am, An);  \
    RUNTIME_ACCESS_WRITE_CACHED(A, Am, An)

#define RUNTIME_RANK_CHANGED(rank)              \
    RUNTIME_PRUNING_STATS_RANK_CHANGED(rank)

#define RUNTIME_END_ACCESS_DECLARATION          \
    RUNTIME_PRUNING_STATS_END_ACCESS_DECLARATION;

#endif /* _HICMA_STARPU_H_ */
