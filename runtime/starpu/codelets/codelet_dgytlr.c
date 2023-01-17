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
#include <hicma.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

DCODELETS_HEADER(gytlr)

#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday
#include <stdlib.h>//FIXME for malloc

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
extern int print_index;
extern int store_only_diagonal_tiles;
extern int global_check;
extern int print_mat;
extern void _printmat(double * A, int m, int n, int ld);

int gytlr_tile_ii = -1;
int gytlr_tile_jj = -1;

void dgytlr( int m, int n, /*dimension of squareAD*/
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
        )
{
    if(gytlr_tile_ii >= 0) {
        ii = gytlr_tile_ii;
        printf("%s %d: Using fixed i:%d\n", __FILE__, __LINE__, ii);
    }
    if(gytlr_tile_jj >= 0) {
        jj = gytlr_tile_jj;
        printf("%s %d: Using fixed j:%d\n", __FILE__, __LINE__, jj);
    }
    int64_t i, j;
    //printf("m:%d n:%d bigM:%d m0:%d n0:%d\n", m, n, bigM, m0, n0);
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    if(print_index){
        fprintf(stderr, "+GYTLR\t|(%d,%d) m:%d n:%d lda:%d ldu:%d ldv:%d\n", ii, jj, m, n, lda, ldu, ldv);
    }

    int shape[2];
    int rank = 0;
    int oversample = 10;
    double *work;
    int *iwork;
    STARSH_blrf* blrf =  HICMA_get_starsh_format();
    STARSH_cluster *RC = blrf->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;
    double *saveAD;
                                                    // allocate space for dense tile
    if((ii != jj && store_only_diagonal_tiles == 1) // if tile is off diagonal and
                                                    // and only diagonal tiles are stored in a tall and skinny matrix
                                                    // store_only_diagonal_tiles is here because
                                                    // code can use AD to store dense tiles
            ||                                      // OR
            compress_diag == 1) {                   // diagonals are also compressed so AD may not used perhaps
        saveAD = AD;
        //AD = malloc(sizeof(double) * m * n);
        AD = malloc(sizeof(double) * lda * n);
        assert(m==lda);
    }

    blrf->problem->kernel(m, n, RC->pivot+RC->start[ii], CC->pivot+CC->start[jj],
            RD, CD, AD, lda);

/*    {*/
        /*if (ii != jj || compress_diag == 1) { */
            /*if(store_only_diagonal_tiles == 1) {*/
                /*assert(AD != saveAD);*/
                /*free(AD);*/
            /*}*/
        /*}*/
        /*return; //TODO*/
    /*}*/
    if(global_check == 1){
       char chall = 'A';
       dlacpy_(&chall, &m, &n, AD, &lda, Dense, &lda);
       //printf("Original problem is copied :%d,%d\n", ii,jj);
    }
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
    if (ii != jj || compress_diag == 1) { // do not try to compress diagonal blocks if it is not enforced
        //AD is m x n. AU and AV are m x maxrank and n x maxrank correspondingly
        //starsh_kernel_drsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        //starsh_dense_dlrrsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        starsh_dense_dlrrsdd(m, n, AD, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);
    
        
    if(0)cblas_dgemm( //for testing purposes 
            CblasColMajor,
            CblasNoTrans, CblasTrans,
            m, n, rank,
            1.0,  AU,  ldu,
            AV,  ldv,
            0.0, Dense, lda);

        if(rank == -1){ //means that tile is dense.
            rank = m;
            fprintf(stderr, "%s %s %d: Dense off-diagonal block (%d,%d). maxrank:%d\n", __FILE__, __func__, __LINE__, ii, jj, maxrank);
            exit(0);
        }
        if(rank == 0) rank = 1;
        Ark[0] = rank;
        if(store_only_diagonal_tiles == 1) {
            assert(AD != saveAD);
            free(AD);
        }
        if(print_mat){
            printf("%d\tgytlr-UV-output\n", __LINE__);
            _printmat(AD, m, n, lda);
            _printmat(AU, m, rank, ldu);
            _printmat(AV, ldv, rank, ldv);
        }
    } else {
        Ark[0] = m;
        if(print_mat){
            printf("%d\tgytlr-DENSE-output\n", __LINE__);
            _printmat(AD, m, m, lda);
        }
    }

    /*printf("m:%d n:%d Tile %d,%d rk:%d lda:%d maxrank:%d oversample:%d tol:%.2e AD:%p AU:%p AV:%p\n", */
            /*m, n, ii, jj, rank, lda, maxrank, oversample, tol, AD, AU, AV);*/
    /*double *tmp = AD;*/
    /*for (i = 0; i < m; ++i) {*/
        /*for (j=0; j<n; ++j ) {*/
            /*printf("%+.2e ",tmp[j*lda+i]);*/
        /*}*/
        /*printf("\n");*/
    /*}*/
    /*printf("\n");*/

    free(work);
    free(iwork);
    if(print_index){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "-GYTLR\t|(%d,%d) rk:%g m:%d n:%d lda:%d ldu:%d ldv:%d\t\t\t\t\tGYTLR: %.4f\n",ii,jj,Ark[0],m, n, lda, ldu, ldv,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
}



/*   HICMA_TASK_dgytlr - Generate a tile for random matrix. */

void HICMA_TASK_dgytlr( const HICMA_option_t *options,
                        int m, int n,
                        const HICMA_desc_t *AUV,
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
    struct starpu_codelet *codelet = &cl_dgytlr;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(AUV, Am, An);
    HICMA_ACCESS_W(Ark, Am, An);
    HICMA_ACCESS_RW(Dense, Am, An);
    HICMA_END_ACCESS_DECLARATION;

    //printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);
        //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE,    &m,                      sizeof(int),
            STARPU_VALUE,    &n,                      sizeof(int),
            STARPU_VALUE,    &nAUV,                      sizeof(int),
            STARPU_W,         RTBLKADDR(AUV, double, Am, An),
            STARPU_W,         RTBLKADDR(Ark, double, Am, An),
            STARPU_RW,         RTBLKADDR(Dense, double, Am, An), // _R must be _W SERIOUSLY. BUT _W STALLS ON SHAHEEN. FIXME
            STARPU_VALUE,  &lda,                      sizeof(int),
            STARPU_VALUE,  &ldu,                      sizeof(int),
            STARPU_VALUE,  &ldv,                      sizeof(int),
            STARPU_VALUE, &bigM,                      sizeof(int),
            STARPU_VALUE,   &Am,                      sizeof(int),
            STARPU_VALUE,   &An,                      sizeof(int),
            STARPU_VALUE, &seed,   sizeof(unsigned long long int),
            STARPU_VALUE,   &maxrank,                      sizeof(int),
            STARPU_VALUE,   &tol,                      sizeof(double),
            STARPU_VALUE,   &compress_diag,            sizeof(int),
            STARPU_PRIORITY,    options->priority,
            STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "dgytlr",
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
    double *AD = NULL;
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
    Ark = (double *)STARPU_MATRIX_GET_PTR(descr[1]);
    Dense = (double *)STARPU_MATRIX_GET_PTR(descr[2]);


    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &lda, &ldu, &ldv, &bigM, &m0, &n0, &seed, &maxrank, &tol, &compress_diag );

    double *AU = AUV;
    int nAU = nAUV/2;
    assert(ldu == ldv);
    size_t nelm_AU = (size_t)ldu * (size_t)nAU;
    double *AV = &(AUV[nelm_AU]);

    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,HICMA_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    dgytlr( m, n,
            AU,
            AV,
            AD,
            Ark,
            lda,
            ldu,
            ldv,
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
CODELETS_CPU(dgytlr, 3, cl_dgytlr_cpu_func)
