/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zhagcm.c
 *
 *  HiCMA HCORE routine for generating dense matrix.
 *  
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
//#include "hcore/include/hcore.h"
#include "morse.h"
#include "hcore_z.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday

#include "hicma.h"
#include "hicma_common.h"

#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"

extern void _printmat(double * A, int m, int n, int ld);

void HCORE_zhagcm( int m, int n, /*dimension of squareAD*/
        double *AU,
        double *AV,
        double *Ark,
        int ldu,
        int ldv,
        int tile_row_index, int tile_column_index,
        int maxrank, double tol)
{
    int ii = tile_row_index; 
    int jj = tile_column_index; 
    int64_t i, j;
    int lda = ldu; //FIXME ASSUMPTION
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    if(HICMA_get_print_index() == 1 ){
        fprintf(stderr, "%d+HAGcM\t|(%d,%d) m:%d n:%d lda:%d ldu:%d ldv:%d\n",MORSE_My_Mpi_Rank(), ii, jj, m, n, lda, ldu, ldv);
    }

    int shape[2];
    int rank = 0;
    int oversample = 10;
    double *work;
    int *iwork;
    STARSH_cluster *RC = HICMA_get_starsh_format()->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;
    double *AD;
    AD = malloc(sizeof(double) * lda * n);
    //starsh_blrf_get_block(HICMA_get_starsh_format(), ii, jj, shape, &AD);

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
    if(HICMA_get_always_fixed_rank() == 1) {
        HICMA_set_fixed_rank(rank);
    }
    if(HICMA_get_print_mat() == 1){
        printf("%d\thagcm-UV-output\n", __LINE__);
        _printmat(AD, m, n, lda);
        _printmat(AU, m, rank, ldu);
        _printmat(AV, ldv, rank, ldv);
    }

    free(work);
    free(iwork);
    if(HICMA_get_print_index_end() == 1 || 
            HICMA_get_print_index() == 1 ){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-HAGcM\t|(%d,%d) rk:%g m:%d n:%d lda:%d ldu:%d ldv:%d\t\t\t\t\tHAGcM: %.4f\n",MORSE_My_Mpi_Rank(),ii,jj,Ark[0],m, n, lda, ldu, ldv,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
}



