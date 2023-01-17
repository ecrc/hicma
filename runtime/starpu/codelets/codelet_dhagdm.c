/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dhagdm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(hagdm)

#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday

#include <hicma.h>
#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef MKL
  #include <mkl.h>
  #include <mkl_lapack.h>
  //#pragma message("MKL is used")
#else
  #ifdef ARMPL
    #include <armpl.h>
  #else
    #include <cblas.h>
  #endif
  #ifdef LAPACKE_UTILS
    #include <lapacke_utils.h>
  #endif
  #include <lapacke.h>
  //#pragma message("MKL is NOT used")
#endif

//#warning "An experimental feature is enabled!!!" 
extern int steal_lrtile;

extern void _printmat(double * A, int m, int n, int ld);

void dhagdm( 
        int nrows_Dense,
        int ncols_Dense,
        double *Dense,
        int ld_Dense,
        int tile_row_index,
        int tile_col_index,
        int A_mt
        )
{
    if(steal_lrtile == 1 && tile_row_index == tile_col_index){
        if(tile_row_index == A_mt-1) { // steal tile above
            tile_row_index -= 1;
        } else { // still tile below
            tile_row_index += 1;
        }
    }
    struct timeval tvalBefore, tvalAfter; 
    gettimeofday (&tvalBefore, NULL);
    STARSH_cluster *RC = HICMA_get_starsh_format()->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;
    HICMA_get_starsh_format()->problem->kernel(nrows_Dense, ncols_Dense, 
            RC->pivot+RC->start[tile_row_index], 
            CC->pivot+CC->start[tile_col_index],
            RD, CD, Dense, ld_Dense);
}
/**
 * HICMA_TASK_dhagdm - Generate dense matrix from a problem determined according to current global setting of HiCMA library
 */
void HICMA_TASK_dhagdm( const HICMA_option_t *options,
                        int nrows_Dense, int ncols_Dense,
                        const HICMA_desc_t *Dense,
                        int ld_Dense,
                        int tile_row_index,
                        int tile_col_index,
                        int A_mt
                        )
{
    struct starpu_codelet *codelet = &cl_dhagdm;
    void (*callback)(void*) = NULL;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(Dense, tile_row_index, tile_col_index);
    HICMA_END_ACCESS_DECLARATION;
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &nrows_Dense,                   sizeof(int),
            STARPU_VALUE,    &ncols_Dense,                   sizeof(int),
            STARPU_W,         RTBLKADDR(Dense, double, tile_row_index, tile_col_index),
            STARPU_VALUE,  &ld_Dense,                        sizeof(int),
            STARPU_VALUE,   &tile_row_index,                 sizeof(int),
            STARPU_VALUE,   &tile_col_index,                 sizeof(int),
            STARPU_VALUE,   &A_mt,                           sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dhagdm",
#endif
            0);
}
/**
 * cl_dhagdm_cpu_func - Generate a tile for random matrix.
 */

#if !defined(CHAMELEON_SIMULATION)
static void cl_dhagdm_cpu_func(void *descr[], void *cl_arg)
{
    int nrows_Dense;
    int ncols_Dense;
    int ld_Dense;
    int tile_row_index;
    int tile_col_index;
    int A_mt;
    int maxrank;
    double *Dense;

    Dense = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &nrows_Dense, &ncols_Dense, &ld_Dense, &tile_row_index, &tile_col_index, &A_mt);
    dhagdm(
            nrows_Dense,
            ncols_Dense,
            Dense,
            ld_Dense,
            tile_row_index,
            tile_col_index,
            A_mt
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dhagdm, 1, cl_dhagdm_cpu_func)

DCODELETS_HEADER(hagdmi)
/**
 * HICMA_TASK_dhagdmi - Generate dense matrix from a problem determined according to current global setting of HiCMA library
 * This function takes indices of tiles of problem.
 */
void HICMA_TASK_dhagdmi( const HICMA_option_t *options,
                        int nrows_Dense, int ncols_Dense,
                        const HICMA_desc_t *Dense,
                        int ld_Dense,
                        int tile_row_index,
                        int tile_col_index,
                        int problem_row_index,
                        int problem_col_index
                        )
{
    struct starpu_codelet *codelet = &cl_dhagdmi;
    void (*callback)(void*) = NULL;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(Dense, tile_row_index, tile_col_index);
    HICMA_END_ACCESS_DECLARATION;
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &nrows_Dense,                   sizeof(int),
            STARPU_VALUE,    &ncols_Dense,                   sizeof(int),
            STARPU_W,         RTBLKADDR(Dense, double, tile_row_index, tile_col_index),
            STARPU_VALUE,  &ld_Dense,                        sizeof(int),
            STARPU_VALUE,   &tile_row_index,                 sizeof(int),
            STARPU_VALUE,   &tile_col_index,                 sizeof(int),
            STARPU_VALUE,   &problem_row_index,                 sizeof(int),
            STARPU_VALUE,   &problem_col_index,                 sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dhagdm",
#endif
            0);
}

/** cl_dhagdm_cpu_func - Generate a tile for random matrix. 
 * This function takes indices of tiles of problem.
 */

#if !defined(CHAMELEON_SIMULATION)
static void cl_dhagdmi_cpu_func(void *descr[], void *cl_arg)
{
    int nrows_Dense;
    int ncols_Dense;
    int ld_Dense;
    int tile_row_index;
    int tile_col_index;
    int maxrank;
    double *Dense;
    int problem_row_index;
    int problem_col_index;

    Dense = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &nrows_Dense, &ncols_Dense, &ld_Dense, &tile_row_index, &tile_col_index, &problem_row_index, &problem_col_index);
    dhagdm(
            nrows_Dense,
            ncols_Dense,
            Dense,
            ld_Dense,
            problem_row_index,
            problem_col_index, -1
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dhagdmi, 1, cl_dhagdmi_cpu_func)
