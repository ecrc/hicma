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

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(hagcm)


#include <assert.h>
#include <stdio.h>
#include <stdlib.h>  //FIXME remove mallocs from this code
#include <sys/time.h>//FIXME for gettimeofday

#include "hicma.h"
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
int steal_lrtile = 0; //non-zero values are for experimental reasons. Otherwise set to 0!!!!!


extern void _printmat(double * A, int m, int n, int ld);

void zhagcm( int m, int n, /*dimension of squareAD*/
        double *AU,
        double *AV,
        double *Ark,
        int ldu,
        int ldv,
        int tile_row_index, int tile_column_index,
        int maxrank, double tol,
        int A_mt)
{
    int ii = tile_row_index; 
    int jj = tile_column_index; 
    if(steal_lrtile == 1 && ii == jj){
        if(ii == A_mt-1) { // steal tile above
            ii = ii - 1;
        } else { // still tile below
            ii = ii + 1;
        }
    }
    int64_t i, j;
    int lda = ldu; //FIXME ASSUMPTION
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);

    int shape[2];
    int rank = 0;
    int oversample = 10;
    double *work;
    int *iwork;
    STARSH_cluster *RC = HICMA_get_starsh_format()->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;
    double *AD;
    AD = malloc(sizeof(double) * lda * n);

    HICMA_get_starsh_format()->problem->kernel(m, n, RC->pivot+RC->start[ii], CC->pivot+CC->start[jj],
            RD, CD, AD, lda);
    int mn = m;
    int mn2 = maxrank+oversample;
    if(mn2 > mn)
        mn2 = mn;
    // Get size of temporary arrays
    size_t lwork = n, lwork_sdd = (4*mn2+7)*mn2;
    if(lwork_sdd > lwork)
        lwork = lwork_sdd;
    lwork += (size_t)mn2*(2*n+m+mn2+1);
    size_t liwork = 8*mn2;
    // Allocate temporary arrays
    //STARSH_MALLOC(iwork, liwork);
    iwork = malloc(sizeof(*iwork) * liwork);
    if(iwork == NULL) {
        fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! liwork:%d", __FILE__, __func__, __LINE__, liwork);
        exit(-1);
    }
    //STARSH_MALLOC(work, lwork);
    work = malloc(sizeof(*work) * lwork);
    if(work == NULL) {
        fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! lwork:%d", __FILE__, __func__, __LINE__, lwork);
        exit(-1);
    }
    starsh_dense_dlrrsdd(m, n, AD, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);

    Ark[0] = rank;

    free(work);
    free(iwork);
}
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
                        int maxrank, double tol,
                        int A_mt
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
            STARPU_VALUE,   &A_mt,                      sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zhagcm",
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
    int A_mt;

    AUV = (double *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[1]);


    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &ldu, &ldv, &tile_row_index, &tile_column_index, &maxrank, &tol, &A_mt);

    double *AU = AUV;
    int nAU = nAUV/2;
    assert(ldu == ldv);
    size_t nelm_AU = (size_t)ldu * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);

    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,MORSE_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    zhagcm( m, n,
            AU,
            AV,
            Ark,
            ldu,
            ldv,
            tile_row_index, tile_column_index,
            maxrank, tol, A_mt
            );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zhagcm, 2, cl_zhagcm_cpu_func)
