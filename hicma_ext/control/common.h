/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file common.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2015 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon common header file
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
/**
 *  HICMA facilities of interest to both HICMA core developer
 *  and also of interest to HICMA community contributor.
 */
#ifndef _HICMA_CHAM_COMMON_H_
#define _HICMA_CHAM_COMMON_H_


#if defined( _WIN32 ) || defined( _WIN64 )
#include <io.h>
#else
#include <unistd.h>
#endif

/**
 * Implementation headers
 */
#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
#include <cuda.h>
#include <cuda_runtime.h>
#include <cuda_runtime_api.h>
#if defined(CHAMELEON_USE_CUBLAS_V2)
#include <cublas.h>
#include <cublas_v2.h>
#else
#include <cublas.h>
#endif
#endif

#if defined(CHAMELEON_USE_OPENCL) && !defined(CHAMELEON_SIMULATION)
#include <OpenCL/cl.h>
#endif

#if defined(HICMA_USE_MPI)
#include <mpi.h>
#endif

/**
 *  Line to avoid conflict with other linear algebra libraries, because, we
 *  don't know why but lapacke provide a wrong interface of lapack in fortran
 */
#ifndef LAPACK_NAME
#define LAPACK_NAME(a, b) lapackef77_##a
#endif

/**
 *  Chameleon header files
 */

#include <hicma.h>
#include "hicma_global.h"
#include "hicma_auxiliary.h"
#include "hicma_context.h"
#include "hicma_descriptor.h"
#include "hicma_async.h"

/**
 *  Global shortcuts
 */
#define HICMA_RANK        hicma_rank(hicma)
#define HICMA_SIZE        hicma->world_size
#define HICMA_GRPSIZE     hicma->group_size
#define HICMA_NB          hicma->nb
#define HICMA_IB          hicma->ib
#define HICMA_NBNBSIZE    hicma->nbnbsize
#define HICMA_IBNBSIZE    hicma->ibnbsize
#define HICMA_SCHEDULING  hicma->scheduling
#define HICMA_RHBLK       hicma->rhblock
#define HICMA_TRANSLATION hicma->translation
#define HICMA_PARALLEL    hicma->parallel_enabled
#define HICMA_PROFILING   hicma->profiling_enabled
#if defined(HICMA_USE_MPI)
#define HICMA_MPI_RANK    hicma->my_mpi_rank
#define HICMA_MPI_SIZE    hicma->mpi_comm_size
#endif

/**
 *  IPT internal define
 */
#define HicmaIPT_NoDep   0
#define HicmaIPT_Panel   1
#define HicmaIPT_All     2

/**
 *  Global array of LAPACK constants
 */
extern char *hicma_lapack_constants[];
#define hicma_lapack_const(hicma_const) hicma_lapack_constants[hicma_const][0]

#ifdef __cplusplus
extern "C" {
#endif


#include <control/hicma_compute_s.h>
#include <control/hicma_compute_d.h>
#include <control/hicma_compute_c.h>
#include <control/hicma_compute_z.h>

/*
void hicma_pdlag2s(HICMA_context_t *hicma);
void hicma_pzlag2c(HICMA_context_t *hicma);
void hicma_pslag2d(HICMA_context_t *hicma);
void hicma_pclag2z(HICMA_context_t *hicma);
*/

#ifdef __cplusplus
}
#endif

#endif
