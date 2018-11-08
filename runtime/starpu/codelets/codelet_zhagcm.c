/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_zhagcm.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
//#include "runtime/starpu/include/runtime_codelet_z.h"
#include "hcore_z.h"

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(hagcm)


//CHAMELEON_CL_CB(zhagcm,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N)
/**
 * HICMA_TASK_zhagcm - Generate compressed matrix from a problem determined according to current global setting of HiCMA library
 */
void HICMA_TASK_zhagcm( const MORSE_option_t *options,
                        int m, int n,
                        const MORSE_desc_t *AUV,
                        const MORSE_desc_t *Ark,
                        int Am, int An, 
                        int ldu,
                        int ldv,
                        int maxrank, double tol
                        )
{
    struct starpu_codelet *codelet = &cl_zhagcm;
    //void (*callback)(void*) = options->profiling ? cl_zhagcm_callback : NULL;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(AUV, Am, An);
    MORSE_ACCESS_W(Ark, Am, An);
    MORSE_END_ACCESS_DECLARATION;

    //printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);
        //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                      sizeof(int),
            STARPU_VALUE,    &n,                      sizeof(int),
            STARPU_VALUE,    &nAUV,                      sizeof(int),
            STARPU_W,         RTBLKADDR(AUV, double, Am, An),
            STARPU_W,         RTBLKADDR(Ark, double, Am, An),
            STARPU_VALUE,  &ldu,                      sizeof(int),
            STARPU_VALUE,  &ldv,                      sizeof(int),
            STARPU_VALUE,   &Am,                      sizeof(int),
            STARPU_VALUE,   &An,                      sizeof(int),
            STARPU_VALUE,   &maxrank,                      sizeof(int),
            STARPU_VALUE,   &tol,                      sizeof(double),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "hcore_zhagcm",
#endif
            0);
}

/*   cl_zhagcm_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zhagcm_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int nAUV;
    double *AUV;
    double *Ark;
    int ldu;
    int ldv;
    int tile_row_index;
    int tile_column_index;
    int maxrank;
    double tol;

    AUV = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[1]);


    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &ldu, &ldv, &tile_row_index, &tile_column_index, &maxrank, &tol);

    double *AU = AUV;
    int nAU = nAUV/2;
    assert(ldu == ldv);
    size_t nelm_AU = (size_t)ldu * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);

    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,MORSE_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    HCORE_zhagcm( m, n,
            AU,
            AV,
            Ark,
            ldu,
            ldv,
            tile_row_index, tile_column_index,
            maxrank, tol
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhagcm, 2, cl_zhagcm_cpu_func)
