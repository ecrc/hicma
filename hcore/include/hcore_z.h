/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_z.h
 *
 *  HiCMA HCORE kernels
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#ifndef _HICMA_HCORE_Z_H_
#define _HICMA_HCORE_Z_H_

#define COMPLEX

#include "morse.h"
#ifdef __cplusplus
extern "C" {
#endif
    /** ****************************************************************************
     *  Declarations of serial kernels - alphabetical order
     **/
    void HCORE_zgemm(MORSE_enum transA, int transB,
            int M, int N,
            double alpha, 
            double *AU, 
            double *AV, 
            double *Ark, 
            int LDA,
            double *BU,
            double *BV,
            double *Brk,
            int LDB,
            double beta,
            double *CU,
            double *CV,
            double *Crk,
            int LDC,
            int rk,
            int maxrk,
            double acc,
            double* work
                );
    void HCORE_zgemm_fast(MORSE_enum transA, int transB,
            int M, int N, 
            double alpha,
            double *AU,
            double *AV,
            double *Ark,
            int LDA,
            double *BU,
            double *BV,
            double *Brk,
            int LDB,
            double beta,
            double *CU,
            double *CV,
            double *Crk,
            int LDC,
            int rk,
            int maxrk,
            double acc,
            double* work
            );
    void HCORE_zgemmbdcd(MORSE_enum transA, MORSE_enum transB,
            int M, int N,
            double alpha, 
            double *AU, 
            double *AV, 
            double *Ark, 
            int LDA,
            double *BD, 
            int LDB,
            double beta, 
            double *CD, 
            int LDC,
            double *work 
            );
    void HCORE_zgytlr(int m, int n,
            double *AU,
            double *AV,
            double *AD,
            double *Ark,
            int lda,
            int ldu,
            int ldv,
            int bigM, int m0, int n0, unsigned long long int seed,
            int maxrank, double tol,
            int compress_diag,
            double *Dense
            );
    void HCORE_zhagdm( 
            int nrows_Dense,
            int ncols_Dense,
            double *Dense,
            int ld_Dense,
            int tile_row_index,
            int tile_col_index
            );
    void HCORE_zhagcm( int m, int n, 
            double *AU,
            double *AV,
            double *Ark,
            int ldu,
            int ldv,
            int tile_row_index, 
            int tile_column_index, 
            int maxrank, double tol);
    void HCORE_zsyrk(MORSE_enum uplo, MORSE_enum trans,
            int N, int K,
            double alpha,
            const double *AU, int LDAU,
            const double *AV, int LDAV,
            double beta,
            double *CD, int LDCD,
            double* work
            );
    void HCORE_zuncompress(MORSE_enum transA, MORSE_enum transB,
            int M, int N,
            double alpha, 
            double *AU, 
            double *Ark, 
            int LDA,
            double *BV,
            double *Brk,
            int LDB,
            double beta,
            double *CD,
            int LDC
            );
#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
