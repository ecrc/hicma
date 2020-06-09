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

    /** ****************************************************************************
     *  Declarations of serial kernels - alphabetical order
     **/
    
#if defined(HICMA_DOUBLE)
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
            int tile_col_index,
            int A_mt
            );
    void HCORE_zhagcm( int m, int n, 
            double *AU,
            double *AV,
            double *Ark,
            int ldu,
            int ldv,
            int tile_row_index, 
            int tile_column_index, 
            int maxrank, double tol, int A_mt);
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


#undef COMPLEX

#endif

#if defined(HICMA_COMPLEX)
void HCORE_zgemm(MORSE_enum transA, int transB,
                 int M, int N,
                 MORSE_Complex64_t alpha,
                 MORSE_Complex64_t *AU,
                 MORSE_Complex64_t *AV,
                 double *Ark,
                 int LDA,
                 MORSE_Complex64_t *BU,
                 MORSE_Complex64_t *BV,
                 double *Brk,
                 int LDB,
                 MORSE_Complex64_t beta,
                 MORSE_Complex64_t *CU,
                 MORSE_Complex64_t *CV,
                 double *Crk,
                 int LDC,
                 int rk,
                 int maxrk,
                 double acc,
                 MORSE_Complex64_t* work
                 );
    void HCORE_dgemm(MORSE_enum transA, int transB,
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
                      MORSE_Complex64_t alpha,
                      MORSE_Complex64_t *AU,
                      MORSE_Complex64_t *AV,
                      MORSE_Complex64_t *Ark,
                      int LDA,
                      MORSE_Complex64_t *BU,
                      MORSE_Complex64_t *BV,
                      MORSE_Complex64_t *Brk,
                      int LDB,
                      MORSE_Complex64_t beta,
                      MORSE_Complex64_t *CU,
                      MORSE_Complex64_t *CV,
                      MORSE_Complex64_t *Crk,
                      int LDC,
                      int rk,
                      int maxrk,
                      double acc,
                      MORSE_Complex64_t* work
                      );
void HCORE_zgemmbdcd(MORSE_enum transA, MORSE_enum transB,
                     int M, int N,
                     MORSE_Complex64_t alpha,
                     MORSE_Complex64_t *AU,
                     MORSE_Complex64_t *AV,
                     double *Ark,
                     int LDA,
                     MORSE_Complex64_t *BD,
                     int LDB,
                     MORSE_Complex64_t beta,
                     MORSE_Complex64_t *CD,
                     int LDC,
                     MORSE_Complex64_t *work
                     );
void HCORE_zgytlr(int m, int n,
                  MORSE_Complex64_t *AU,
                  MORSE_Complex64_t *AV,
                  MORSE_Complex64_t *AD,
                  double *Ark,
                  int lda,
                  int ldu,
                  int ldv,
                  int bigM, int m0, int n0, unsigned long long int seed,
                  int maxrank, double tol,
                  int compress_diag,
                  MORSE_Complex64_t *Dense
                  );
/*void HCORE_zhagdm(
                  int nrows_Dense,
                  int ncols_Dense,
                  MORSE_Complex64_t *Dense,
                  int ld_Dense,
                  int tile_row_index,
                  int tile_col_index,
                  int A_mt
                  );
void HCORE_zhagcm( int m, int n,
                  MORSE_Complex64_t *AU,
                  MORSE_Complex64_t *AV,
                  MORSE_Complex64_t *Ark,
                  int ldu,
                  int ldv,
                  int tile_row_index,
                  int tile_column_index,
                  int maxrank, double tol, int A_mt);
void HCORE_zsyrk(MORSE_enum uplo, MORSE_enum trans,
                 int N, int K,
                 MORSE_Complex64_t alpha,
                 const MORSE_Complex64_t *AU, int LDAU,
                 const MORSE_Complex64_t *AV, int LDAV,
                 MORSE_Complex64_t beta,
                 MORSE_Complex64_t *CD, int LDCD,
                 MORSE_Complex64_t* work
                 );
*/
void HCORE_zuncompress(MORSE_enum transA, MORSE_enum transB,
                       int M, int N,
                       MORSE_Complex64_t alpha,
                       MORSE_Complex64_t *AU,
                       double *Ark,
                       int LDA,
                       MORSE_Complex64_t *BV,
                       double *Brk,
                       int LDB,
                       MORSE_Complex64_t beta,
                       MORSE_Complex64_t *CD,
                       int LDC
                       );

#endif
#endif
