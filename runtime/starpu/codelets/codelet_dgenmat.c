/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dgenmat.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Rabab Alomairy
 * @date 2020-05-20
 * @precisions normal z -> c d s
 **/

#include <hicma.h>
#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>
#include <starsh.h>

extern int store_only_diagonal_tiles;

DCODELETS_HEADER(genmat)

void HICMA_TASK_dgenmat( const HICMA_option_t *options,
                         HICMA_desc_t *A, int lda, int Am, int An, int m, int n)
{
	struct starpu_codelet *codelet = &cl_dgenmat;
	void (*callback)(void*) = NULL;

	HICMA_BEGIN_ACCESS_DECLARATION;
	HICMA_ACCESS_W(A, Am, An);
	HICMA_END_ACCESS_DECLARATION;

	// printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);

	//printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
	starpu_insert_task(
			starpu_mpi_codelet(codelet),
			STARPU_W,        RTBLKADDR(A, double, Am, An),
			STARPU_VALUE,    &Am,                      sizeof(int),
			STARPU_VALUE,    &An,                      sizeof(int),
                        STARPU_VALUE,    &lda,                      sizeof(int),
                        STARPU_VALUE,    &m,                      sizeof(int),
                        STARPU_VALUE,    &n,                      sizeof(int),
			STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
			STARPU_NAME, "zgenmat",
#endif
			0);
}

/*   cl_dper_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_dgenmat_cpu_func(void *descr[], void *cl_arg)
{
	int i;
	int j;
        int m, n;
	double *A;
        int lda;

	A = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
        
	starpu_codelet_unpack_args(cl_arg, &i, &j, &lda, &m, &n);

       int shape[2];
    int oversample = 10;
    double *work;
    int *iwork;
    STARSH_blrf* blrf =  HICMA_get_starsh_format();
    STARSH_cluster *RC = blrf->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;

    blrf->problem->kernel(m, n, RC->pivot+RC->start[i], CC->pivot+CC->start[j],
            RD, CD, A, lda);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dgenmat, 1, cl_dgenmat_cpu_func)
