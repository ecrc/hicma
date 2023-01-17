/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_dgytlr.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/

#include <runtime/starpu/hicma_starpu.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(gytlrdiag)

extern void dgytlr( int m, int n, /*dimension of squareAD*/
        double *AU,
        double *AV,
        double *AD,
        double *Ark,
        int lda,
        int ldu,
        int ldv,
        int bigM, int ii, int jj, unsigned long long int seed,
        int maxrank, double tol, int compress_diag,
        double *Dense
        );
/*   HICMA_TASK_dgytlr - Generate a tile for random matrix. */

void HICMA_TASK_dgytlr_diag( const HICMA_option_t *options,
                        int m, int n,
                        const HICMA_desc_t *AUV,
                        const HICMA_desc_t *AD, int ADm, int ADn,
                        const HICMA_desc_t *Ark,
                        int Am, int An, 
                        int lda,
                        int ldu,
                        int ldv,
                        int bigM, int m0, int n0, unsigned long long int seed,
                        int maxrank, double tol,
                        int compress_diag,
                        HICMA_desc_t *Dense
                        )
{
    struct starpu_codelet *codelet = &cl_dgytlrdiag;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(AUV, Am, An);
    HICMA_ACCESS_W(AD, ADm, ADn);
    HICMA_ACCESS_W(Ark, Am, An);
    HICMA_ACCESS_RW(Dense, Am, An);
    HICMA_END_ACCESS_DECLARATION;

    // printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);

        //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                      sizeof(int),
            STARPU_VALUE,    &n,                      sizeof(int),
            STARPU_VALUE,    &nAUV,                      sizeof(int),
            STARPU_W,         RTBLKADDR(AUV, double, Am, An),
            STARPU_W,         RTBLKADDR(AD, double, ADm, ADn),
            STARPU_W,         RTBLKADDR(Ark, double, Am, An),
            STARPU_RW,         RTBLKADDR(Dense, double, Am, An), // _R must be _W SERIOUSLY. BUT _W STALLS ON SHAHEEN. FIXME
            STARPU_VALUE,  &lda,                      sizeof(int),
            STARPU_VALUE,  &ldu,                      sizeof(int),
            STARPU_VALUE,  &ldv,                      sizeof(int),
            STARPU_VALUE, &bigM,                      sizeof(int),
            STARPU_VALUE,   &Am,                      sizeof(int),
            STARPU_VALUE,   &An,                      sizeof(int),
            STARPU_VALUE, &seed,   sizeof(unsigned long long int),
            STARPU_VALUE,   &maxrank,                  sizeof(int),
            STARPU_VALUE,   &tol,                      sizeof(double),
            STARPU_VALUE,   &compress_diag,            sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dgytlr_diag",
#endif
            0);
}

/*   cl_dgytlr_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_dgytlr_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int nAUV;
    double *AUV;
    double *AD;
    double *Ark;
    double *Dense;
    int lda;
    int ldu;
    int ldv;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;
    int maxrank;
    double tol;
    int compress_diag;

    AUV = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    AD = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[2]);
    Dense = (double *)STARPU_MATRIX_GET_PTR(descr[3]);


    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &lda, &ldu, &ldv, &bigM, &m0, &n0, &seed, &maxrank, &tol, &compress_diag );

    double *AU = AUV;
    int nAU = nAUV/2;
    size_t nelm_AU = (size_t)lda * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);

    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,HICMA_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    dgytlr( m, n,
            AU,
            AV,
            AD,
            Ark,
            lda, ldu, ldv, 
            bigM, m0, n0, seed,
            maxrank, tol,
            compress_diag,
            Dense
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(dgytlrdiag, 4, cl_dgytlr_cpu_func)
