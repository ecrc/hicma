/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 * @file hicma_runtime_z.h
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
#ifndef _RUNTIME_ZHCORE_H_
#define _RUNTIME_ZHCORE_H_

#ifdef __cplusplus
extern "C" {
#endif


/** ****************************************************************************
 *  Declarations of HCORE wrappers (called by HiCMA) - alphabetical order
 **/
void HICMA_TASK_hcore_zgemm(const HICMA_option_t *options,
                            HICMA_enum transA, int transB,
                            int m, int n,
                            HICMA_Complex64_t alpha,
                            const HICMA_desc_t *AUV,
                            const HICMA_desc_t *Ark,
                            int Am, int An, int lda,
                            const HICMA_desc_t *BUV,
                            const HICMA_desc_t *Brk,
                            int Bm, int Bn, int ldb,
                            HICMA_Complex64_t beta,
                            const HICMA_desc_t *CUV,
                            const HICMA_desc_t *Crk,
                            int Cm, int Cn, int ldc,
                            int rk, int maxrk,
                            double acc);

void HICMA_TASK_zgemm_bdcd(const HICMA_option_t *options,
                           HICMA_enum transA, int transB,
                           int m, int n,
                           double alpha,
                           const HICMA_desc_t *AUV,
                           const HICMA_desc_t *Ark,
                           int Am, int An, int lda,
                           const HICMA_desc_t *BD,
                           int Bm, int Bn, int ldb,
                           double beta,
                           const HICMA_desc_t *CD,
                           int Cm, int Cn, int ldc
);

void HICMA_TASK_zgytlr_diag(const HICMA_option_t *options,
                            int m, int n,
                            const HICMA_desc_t *AUV,
                            const HICMA_desc_t *AD, int ADm, int ADn,
                            const HICMA_desc_t *Ark,
                            int Am, int An,
                            int lda, int ldu, int ldv,
                            int bigM, int m0, int n0, unsigned long long int seed,
                            int maxrank, double tol, int compress_diag,
                            const HICMA_desc_t *Dense
);

void HICMA_TASK_zgytlr(const HICMA_option_t *options,
                       int m, int n,
                       const HICMA_desc_t *AUV,
                       const HICMA_desc_t *Ark,
                       int Am, int An,
                       int lda, int ldu, int ldv,
                       int bigM, int m0, int n0, unsigned long long int seed,
                       int maxrank, double tol, int compress_diag,
                       const HICMA_desc_t *Dense
);

void HICMA_TASK_zhagdm(const HICMA_option_t *options,
                       int nrows_Dense, int ncols_Dense,
                       const HICMA_desc_t *Dense,
                       int ld_Dense,
                       int tile_row_index,
                       int tile_col_index,
                       int A_mt
);

void HICMA_TASK_zhagdmi(const HICMA_option_t *options,
                        int nrows_Dense, int ncols_Dense,
                        const HICMA_desc_t *Dense,
                        int ld_Dense,
                        int tile_row_index,
                        int tile_col_index,
                        int problem_row_index,
                        int problem_col_index
);

void HICMA_TASK_zhagcm(const HICMA_option_t *options,
                       int m, int n,
                       const HICMA_desc_t *AUV,
                       const HICMA_desc_t *Ark,
                       int Am, int An,
                       int ldu,
                       int ldv,
                       int maxrank, double tol,
                       int A_mt
);

void HICMA_TASK_zpotrf(const HICMA_option_t *options,
                       HICMA_enum uplo, int n, int nb,
                       const HICMA_desc_t *A, int Am, int An, int lda,
                       int iinfo);

void HICMA_TASK_zsyrk(const HICMA_option_t *options,
                      HICMA_enum uplo, HICMA_enum trans,
                      int n, int nb,
                      HICMA_Complex64_t alpha,
                      const HICMA_desc_t *AUV, int ldauv,
                      const HICMA_desc_t *Ark,
                      int Am, int An,
                      double beta,
                      const HICMA_desc_t *CD, int ldcd,
                      int Cm, int Cn);

void HICMA_TASK_zgemm_cd(const HICMA_option_t *options,
                         int n, int nb,
                         HICMA_Complex64_t alpha,
                         const HICMA_desc_t *AUV, int ldauv,
                         const HICMA_desc_t *Ark,
                         int Am, int An,
                         const HICMA_desc_t *BUV, int ldbuv,
                         const HICMA_desc_t *Brk,
                         int Bm, int Bn,
                         double beta,
                         const HICMA_desc_t *CD, int ldcd,
                         int Cm, int Cn);

void
HICMA_TASK_hcore_ztrsm(const HICMA_option_t *options, HICMA_enum side, HICMA_enum uplo, HICMA_enum transA, HICMA_enum diag,
                 int m, HICMA_Complex64_t alpha, const HICMA_desc_t *A, int Am, int An, int lda,
                 const HICMA_desc_t *BUV, int Bm, int Bn, int ldb, const HICMA_desc_t *Brk);

void HICMA_TASK_zuncompress(const HICMA_option_t *options,
                            HICMA_enum transA, int transB,
                            int m, int n,
                            HICMA_Complex64_t alpha,
                            const HICMA_desc_t *AUBV,
                            const HICMA_desc_t *Ark,
                            int Am, int An, int lda,
                            double beta,
                            const HICMA_desc_t *CD,
                            int Cm, int Cn, int ldc);

void HICMA_TASK_zlacpyx(const HICMA_option_t *options,
                        HICMA_enum uplo, int m, int n, int nb,
                        int displA, const HICMA_desc_t *A, int Am, int An, int lda,
                        int displB, const HICMA_desc_t *B, int Bm, int Bn, int ldb);

void HICMA_TASK_zlacpy(const HICMA_option_t *options,
                             HICMA_enum uplo, int m, int n, int nb,
                             const HICMA_desc_t *A, int Am, int An, int lda,
                             const HICMA_desc_t *B, int Bm, int Bn, int ldb);

void HICMA_TASK_zlaset(const HICMA_option_t *options,
                       HICMA_enum uplo, int M, int N,
                       HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                       const HICMA_desc_t *A, int Am, int An, int LDA);

void HICMA_TASK_zplrnt(const HICMA_option_t *options,
                       int m, int n, const HICMA_desc_t *A, int Am, int An, int lda,
                       int bigM, int m0, int n0, unsigned long long int seed);

void HICMA_TASK_ztrsm(const HICMA_option_t *options, HICMA_enum side, HICMA_enum uplo, HICMA_enum transA,
                            HICMA_enum diag, int m, int n, HICMA_Complex64_t alpha, const HICMA_desc_t *A, int Am,
                            int An, int lda, const HICMA_desc_t *B, int Bm, int Bn, int ldb);

void HICMA_TASK_zgemm(const HICMA_option_t *options,
                      HICMA_enum transA, int transB,
                      int m, int n, int k,
                      HICMA_Complex64_t alpha,
                      const HICMA_desc_t *A,
                      int Am, int An, int lda,
                      const HICMA_desc_t *B,
                      int Bm, int Bn, int ldb,
                      HICMA_Complex64_t beta,
                      const HICMA_desc_t *C,
                      int Cm, int Cn, int ldc);

void HICMA_TASK_zgenrhs(const HICMA_option_t *options,
                        int m, int n,
                        const HICMA_desc_t *A, int Am, int An,
                        int lda,
                        int bigM, int m0, int n0);

void HICMA_TASK_zgetrf(const HICMA_option_t *options, int n, int nb,
                       const HICMA_desc_t *A, int Am, int An, int lda,
                       int iinfo);

void HICMA_TASK_zgenmat(const HICMA_option_t *options,
                        HICMA_desc_t *A, int lda, int Am, int An, int m, int n);

void
HICMA_TASK_hcore_ztrsmu(const HICMA_option_t *options, HICMA_enum side, HICMA_enum uplo, HICMA_enum transA, HICMA_enum diag,
                  int m, HICMA_Complex64_t alpha, const HICMA_desc_t *A, int Am, int An, int lda,
                  const HICMA_desc_t *BUV, int Bm, int Bn, int ldb, const HICMA_desc_t *Brk);

#ifdef __cplusplus
}
#endif

#endif
