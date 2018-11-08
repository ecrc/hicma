/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_zhagdm.c
 *
 * Codelet for generating dense matrix from a problem determined according to current global setting of HiCMA library.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
#include "hcore_z.h"

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(hagdm)

/**
 * HICMA_TASK_zhagdm - Generate dense matrix from a problem determined according to current global setting of HiCMA library
 */
void HICMA_TASK_zhagdm( const MORSE_option_t *options,
                        int nrows_Dense, int ncols_Dense,
                        const MORSE_desc_t *Dense, 
                        int ld_Dense,
                        int tile_row_index,
                        int tile_col_index
                        )
{
    struct starpu_codelet *codelet = &cl_zhagdm;
    void (*callback)(void*) = NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(Dense, tile_row_index, tile_col_index);
    MORSE_END_ACCESS_DECLARATION;
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &nrows_Dense,                   sizeof(int),
            STARPU_VALUE,    &ncols_Dense,                   sizeof(int),
            STARPU_W,         RTBLKADDR(Dense, double, tile_row_index, tile_col_index),
            STARPU_VALUE,  &ld_Dense,                        sizeof(int),
            STARPU_VALUE,   &tile_row_index,                 sizeof(int),
            STARPU_VALUE,   &tile_col_index,                 sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_zhagdm",
#endif
            0);
}
/**
 * cl_zhagdm_cpu_func - Generate a tile for random matrix.
 */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhagdm_cpu_func(void *descr[], void *cl_arg)
{
    int nrows_Dense;
    int ncols_Dense;
    int ld_Dense;
    int tile_row_index;
    int tile_col_index;
    int maxrank;
    double *Dense;

    Dense = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &nrows_Dense, &ncols_Dense, &ld_Dense, &tile_row_index, &tile_col_index);
    HCORE_zhagdm(
            nrows_Dense,
            ncols_Dense,
            Dense,
            ld_Dense,
            tile_row_index,
            tile_col_index
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhagdm, 1, cl_zhagdm_cpu_func)

ZCODELETS_HEADER(hagdmi)
/**
 * HICMA_TASK_zhagdmi - Generate dense matrix from a problem determined according to current global setting of HiCMA library
 * This function takes indices of tiles of problem.
 */
void HICMA_TASK_zhagdmi( const MORSE_option_t *options,
                        int nrows_Dense, int ncols_Dense,
                        const MORSE_desc_t *Dense, 
                        int ld_Dense,
                        int tile_row_index,
                        int tile_col_index,
                        int problem_row_index,
                        int problem_col_index
                        )
{
    struct starpu_codelet *codelet = &cl_zhagdmi;
    void (*callback)(void*) = NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(Dense, tile_row_index, tile_col_index);
    MORSE_END_ACCESS_DECLARATION;
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
            STARPU_NAME, "hcore_zhagdm",
#endif
            0);
}

/** cl_zhagdm_cpu_func - Generate a tile for random matrix. 
 * This function takes indices of tiles of problem.
 */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhagdmi_cpu_func(void *descr[], void *cl_arg)
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
    HCORE_zhagdm(
            nrows_Dense,
            ncols_Dense,
            Dense,
            ld_Dense,
            problem_row_index,
            problem_col_index
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhagdmi, 1, cl_zhagdmi_cpu_func)
