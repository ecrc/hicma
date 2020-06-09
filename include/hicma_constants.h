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
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 */

#ifndef _HICMA_CONSTANTS_H_
#define _HICMA_CONSTANTS_H_




#define HICMA_STARSH_PROB_RND    1
#define HICMA_STARSH_PROB_SS     2
#define HICMA_STARSH_PROB_RNDUSR 3
#define HICMA_STARSH_PROB_FILE   4
#define HICMA_STARSH_PROB_GEOSTAT   5
#define HICMA_STARSH_PROB_EDSIN  6
#define HICMA_STARSH_PROB_GEOSTAT_POINT   7
#define HICMA_STARSH_PROB_ST_3D_EXP  8
#define HICMA_STARSH_PROB_ST_3D_SQEXP  9

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


char strmatfile[LEN_STR_MAT_FILE];

#endif
