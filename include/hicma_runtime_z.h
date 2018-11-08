/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hicma_runtime_z.h
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
#ifndef _RUNTIME_ZHCORE_H_
#define _RUNTIME_ZHCORE_H_

#define COMPLEX

#ifdef __cplusplus
extern "C" {
#endif


    /** ****************************************************************************
     *  Declarations of HCORE wrappers (called by HiCMA) - alphabetical order
     **/
    void HICMA_TASK_zgemm(const MORSE_option_t *options,
            MORSE_enum transA, int transB,
            int m, int n,
            double alpha, 
            const MORSE_desc_t *AUV, 
            const MORSE_desc_t *Ark, 
            int Am, int An, int lda,
            const MORSE_desc_t *BUV, 
            const MORSE_desc_t *Brk, 
            int Bm, int Bn, int ldb,
            double beta,  
            const MORSE_desc_t *CUV,
            const MORSE_desc_t *Crk, 
            int Cm, int Cn, int ldc,
            int rk, int maxrk,
            double acc);
    void HICMA_TASK_zgemm_bdcd(const MORSE_option_t *options,
            MORSE_enum transA, int transB,
            int m, int n,
            double alpha,
            const MORSE_desc_t *AUV,
            const MORSE_desc_t *Ark,
            int Am, int An, int lda,
            const MORSE_desc_t *BD,
            int Bm, int Bn, int ldb,
            double beta,
            const MORSE_desc_t *CD,
            int Cm, int Cn, int ldc
            );
    void HICMA_TASK_zgytlr_diag( const MORSE_option_t *options,
            int m, int n, 
            const MORSE_desc_t *AUV, 
            const MORSE_desc_t *AD, int ADm, int ADn, 
            const MORSE_desc_t *Ark, 
            int Am, int An, 
            int lda, int ldu, int ldv,
            int bigM, int m0, int n0, unsigned long long int seed,
            int maxrank, double tol, int compress_diag,
            const MORSE_desc_t *Dense 
            );
    void HICMA_TASK_zgytlr( const MORSE_option_t *options,
            int m, int n, 
            const MORSE_desc_t *AUV, 
            const MORSE_desc_t *Ark, 
            int Am, int An, 
            int lda, int ldu, int ldv,
            int bigM, int m0, int n0, unsigned long long int seed,
            int maxrank, double tol, int compress_diag,
            const MORSE_desc_t *Dense 
            );
    void HICMA_TASK_zhagdm( const MORSE_option_t *options,
            int nrows_Dense, int ncols_Dense,
            const MORSE_desc_t *Dense, 
            int ld_Dense,
            int tile_row_index,
            int tile_col_index
            );
    void HICMA_TASK_zhagdmi( const MORSE_option_t *options,
            int nrows_Dense, int ncols_Dense,
            const MORSE_desc_t *Dense, 
            int ld_Dense,
            int tile_row_index,
            int tile_col_index,
            int problem_row_index,
            int problem_col_index
            );
    void HICMA_TASK_zhagcm( const MORSE_option_t *options,
            int m, int n,
            const MORSE_desc_t *AUV,
            const MORSE_desc_t *Ark,
            int Am, int An, 
            int ldu,
            int ldv,
            int maxrank, double tol
            );
    void HICMA_TASK_zpotrf(const MORSE_option_t *options,
            MORSE_enum uplo, int n, int nb,
            const MORSE_desc_t *A, int Am, int An, int lda,
            int iinfo);
    void HICMA_TASK_zsyrk(const MORSE_option_t *options,
            MORSE_enum uplo, MORSE_enum trans,
            int n, int nb,
            double alpha, 
            const MORSE_desc_t *AUV, int ldauv,
            const MORSE_desc_t *Ark,
            int Am, int An, 
            double beta, 
            const MORSE_desc_t *CD, int ldcd,
            int Cm, int Cn);
    void HICMA_TASK_ztrsm(const MORSE_option_t *options,
            MORSE_enum side, MORSE_enum uplo, MORSE_enum transA, MORSE_enum diag,
            int m,
            double alpha, const MORSE_desc_t *A, int Am, int An, int lda,
            const MORSE_desc_t *BUV, int Bm, int Bn, int ldb, const MORSE_desc_t *Brk);
    void HICMA_TASK_zuncompress(const MORSE_option_t *options,
            MORSE_enum transA, int transB,
            int m, int n,
            double alpha, 
            const MORSE_desc_t *AUBV, 
            const MORSE_desc_t *Ark, 
            int Am, int An, int lda,
            double beta,  
            const MORSE_desc_t *CD, 
            int Cm, int Cn, int ldc);
#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
