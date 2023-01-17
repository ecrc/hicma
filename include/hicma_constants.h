/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/*
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/*
 * @file hicma_constants.c
 *
 *  HiCMA constants
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 */

#ifndef _HICMA_CONSTANTS_H_
#define _HICMA_CONSTANTS_H_


/**
 *
 * @brief Chameleon constants - CBLAS & LAPACK
 *  The naming and numbering is consistent with:
 *
 *    1) CBLAS from Netlib (http://www.netlib.org/blas/blast-forum/cblas.tgz),
 *    2) C Interface to LAPACK from Netlib (http://www.netlib.org/lapack/lapwrapc/).
 *
 */
#define HicmaByte              0
#define HicmaInteger           1
#define HicmaRealFloat         2
#define HicmaRealDouble        3
#define HicmaComplexFloat      4
#define HicmaComplexDouble     5

#define HicmaCM              101
#define HicmaRM              102
#define HicmaCCRB            103
#define HicmaCRRB            104
#define HicmaRCRB            105
#define HicmaRRRB            106

#define HicmaNoTrans         111
#define HicmaTrans           112
#define HicmaConjTrans       113

#define HicmaUpper           121
#define HicmaLower           122
#define HicmaUpperLower      123

#define HicmaNonUnit         131
#define HicmaUnit            132

#define HicmaLeft            141
#define HicmaRight           142

#define HicmaOneNorm         171
#define HicmaRealOneNorm     172
#define HicmaTwoNorm         173
#define HicmaFrobeniusNorm   174
#define HicmaInfNorm         175
#define HicmaRealInfNorm     176
#define HicmaMaxNorm         177
#define HicmaRealMaxNorm     178

#define HicmaDistUniform     201
#define HicmaDistSymmetric   202
#define HicmaDistNormal      203

#define HicmaHermGeev        241
#define HicmaHermPoev        242
#define HicmaNonsymPosv      243
#define HicmaSymPosv         244

#define HicmaNoPacking       291
#define HicmaPackSubdiag     292
#define HicmaPackSupdiag     293
#define HicmaPackColumn      294
#define HicmaPackRow         295
#define HicmaPackLowerBand   296
#define HicmaPackUpeprBand   297
#define HicmaPackAll         298

#define HicmaNoVec           301
#define HicmaVec             302
#define HicmaIvec            303

#define HicmaForward         391
#define HicmaBackward        392

#define HicmaColumnwise      401
#define HicmaRowwise         402
#define HicmaTrd            1001
#define HicmaBrd            1002

#define HicmaW               501
#define HicmaA2              502

#define hicma_const_neg(const) (((const-1)^0x01)+1)

/**
 *  HICMA constants - boolean
 */
#define HICMA_FALSE  0
#define HICMA_TRUE   1

#define HICMA_CPU    ((1ULL)<<1)
#define HICMA_CUDA   ((1ULL)<<3)

/**
 *  State machine switches
 */
#define HICMA_WARNINGS        1
#define HICMA_ERRORS          2
#define HICMA_AUTOTUNING      3
#define HICMA_DAG             4
#define HICMA_PROFILING_MODE  5
#define HICMA_PARALLEL_MODE   6
#define HICMA_BOUND           7
#define HICMA_PROGRESS        8
#define HICMA_GEMM3M          9

/**
 *  HICMA constants - configuration parameters
 */
#define HICMA_CONCURRENCY       1
#define HICMA_TILE_SIZE         2
#define HICMA_INNER_BLOCK_SIZE  3
#define HICMA_HOUSEHOLDER_MODE  5
#define HICMA_HOUSEHOLDER_SIZE  6
#define HICMA_TRANSLATION_MODE  7

#define HICMA_FLAT_HOUSEHOLDER  1
#define HICMA_TREE_HOUSEHOLDER  2

#define HICMA_INPLACE           1
#define HICMA_OUTOFPLACE        2

/**
 *  HICMA constants - success & error codes
 */
#define HICMA_SUCCESS                 0
#define HICMA_ERR_NOT_INITIALIZED  -101
#define HICMA_ERR_REINITIALIZED    -102
#define HICMA_ERR_NOT_SUPPORTED    -103
#define HICMA_ERR_ILLEGAL_VALUE    -104
#define HICMA_ERR_NOT_FOUND        -105
#define HICMA_ERR_OUT_OF_RESOURCES -106
#define HICMA_ERR_INTERNAL_LIMIT   -107
#define HICMA_ERR_UNALLOCATED      -108
#define HICMA_ERR_FILESYSTEM       -109
#define HICMA_ERR_UNEXPECTED       -110
#define HICMA_ERR_SEQUENCE_FLUSHED -111

/**
 * Kernels options
 */
#define HICMA_PRIORITY_MIN  0
#define HICMA_PRIORITY_MAX  INT_MAX


/**
 *  Scheduler properties
 */
#define PRIORITY        16
#define CALLBACK        17
#define REDUX           18

/**
 *  HICMA ???
 */
#define HICMA_REQUEST_INITIALIZER {HICMA_SUCCESS}

#define HICMA_STARSH_PROB_RND    1
#define HICMA_STARSH_PROB_SS     2
#define HICMA_STARSH_PROB_RNDUSR 3
#define HICMA_STARSH_PROB_FILE   4
#define HICMA_STARSH_PROB_GEOSTAT   5
#define HICMA_STARSH_PROB_EDSIN  6
#define HICMA_STARSH_PROB_GEOSTAT_POINT   7
#define HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE 8
#define HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT 9
#define HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE 10
#define HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT 11
#define HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN   12
#define HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN_POINT   13
/*
TODO use enums 
//! Enum for backend types
enum STARSH_BACKEND
{
    STARSH_BACKEND_NOTSELECTED = -2,
    //!< Backend has not been yet selected
    STARSH_BACKEND_NOTSUPPORTED = -1,
    //!< Error, backend is not supported
    STARSH_BACKEND_SEQUENTIAL = 0,
    //!< Sequential
    STARSH_BACKEND_OPENMP = 1,
    //!< OpenMP
    STARSH_BACKEND_MPI = 2,
    //!< MPI
    STARSH_BACKEND_MPI_OPENMP = 3,
    //!< Hybrid MPI + OpenMP
    STARSH_BACKEND_STARPU = 4,
    //!< StarPU (without MPI)
    STARSH_BACKEND_MPI_STARPU = 5
    //!< StarPU (with MPI)
};
 */

#define LEN_STR_MAT_FILE       512


extern char strmatfile[LEN_STR_MAT_FILE];

#endif
