/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dgenrhs.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Rabab Alomairy
 * @date 2020-05-20
 * @precisions normal z -> c d s
 **/

#include <starsh-rbf.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(genrhs)

void HICMA_TASK_dgenrhs( const HICMA_option_t *options,
                        int m, int n,
                        const HICMA_desc_t *A, int Am, int An,
                        int lda,
                        int bigM, int m0, int n0
                        )
{
    struct starpu_codelet *codelet = &cl_dgenrhs;
  
    void (*callback)(void*) = NULL;
    int nb = A->nb;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(A, Am, An);
    HICMA_END_ACCESS_DECLARATION;


    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                      sizeof(int),
            STARPU_VALUE,    &n,                      sizeof(int),
            STARPU_W,         RTBLKADDR(A, double, Am, An),
            STARPU_VALUE,  &lda,                      sizeof(int),
            STARPU_VALUE, &bigM,                      sizeof(int),
            STARPU_VALUE,   &m0,                      sizeof(int),
            STARPU_VALUE,   &n0,                      sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zgenrhs",
#endif
            0);
}

/*   cl_dgenrhs_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_dgenrhs_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    double *A, *mesh;
    int lda;
    int bigM;
    int m0;
    int n0;
    A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);

    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &bigM, &m0, &n0);


     starsh_generate_3d_virus_rhs( m, A);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dgenrhs, 1, cl_dgenrhs_cpu_func)
