/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zgytlr.c
 *
 *  HiCMA HCORE routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
//#include "hcore/include/hcore.h"
#include "morse.h"
#include "hcore_z.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday
#include <stdlib.h>//FIXME for malloc

#define COMPLEX
#undef REAL
//  HCORE_zgytlr - Generate a tile for random matrix.

#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"
extern STARSH_blrf *mpiF;
extern int print_index;
int print_index_end;
extern int store_only_diagonal_tiles;
extern int global_check;
extern int print_mat;
extern void _printmat(double * A, int m, int n, int ld);

//extern int global_always_fixed_rank; //FIXME this does not work, I do not know why
//extern int global_fixed_rank;
int global_always_fixed_rank;
int global_fixed_rank;
int gytlr_tile_ii = -1;
int gytlr_tile_jj = -1;

void HCORE_zgytlr( int m, int n, /*dimension of squareAD*/
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
        fprintf(stderr, "%d+GYTLR\t|(%d,%d) m:%d n:%d lda:%d ldu:%d ldv:%d\n",MORSE_My_Mpi_Rank(), ii, jj, m, n, lda, ldu, ldv);
    }

    int shape[2];
    int rank = 0;
    int oversample = 10;
    double *work;
    int *iwork;
    STARSH_cluster *RC = mpiF->row_cluster, *CC = RC;
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
    //starsh_blrf_get_block(mpiF, ii, jj, shape, &AD);

    mpiF->problem->kernel(m, n, RC->pivot+RC->start[ii], CC->pivot+CC->start[jj],
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
       LAPACK_dlacpy(&chall, &m, &n, AD, &lda, Dense, &lda);
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
            fprintf(stderr, "%s %s %d: Dense off-diagonal block (%d,%d)\n", __FILE__, __func__, __LINE__, ii, jj);
            exit(0);
        }
        Ark[0] = rank;
        if(store_only_diagonal_tiles == 1) {
            assert(AD != saveAD);
            free(AD);
        }
        if(global_always_fixed_rank == 1) {
            global_fixed_rank = rank;
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
    if(print_index || print_index_end ){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-GYTLR\t|(%d,%d) rk:%g m:%d n:%d lda:%d ldu:%d ldv:%d\t\t\t\t\tGYTLR: %.4f\n",MORSE_My_Mpi_Rank(),ii,jj,Ark[0],m, n, lda, ldu, ldv,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
}



