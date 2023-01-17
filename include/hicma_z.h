/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file hicma_z.h
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/

#ifndef _HICMA_Z_H_
#define _HICMA_Z_H_

#undef REAL
#define COMPLEX

#include <hicma.h>

#ifdef __cplusplus
extern "C" {
#endif


//FIXME Naive interfaces taking only arrays are not implemented yet
int HICMA_zpotrf(HICMA_enum uplo, int N, double *A, int LDA);

int HICMA_zpotrf_Tile(HICMA_enum uplo,
                      HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark,
                      int rk, int maxrk, double acc);

int HICMA_zpotrf_Tile_Async(HICMA_enum uplo,
                            HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark,
                            int rk, int maxrk, double acc,
                            HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zgemm_Tile(HICMA_enum transA, HICMA_enum transB,
                     double alpha,
                     HICMA_desc_t *AUV, HICMA_desc_t *Ark,
                     HICMA_desc_t *BUV, HICMA_desc_t *Brk,
                     double beta,
                     HICMA_desc_t *CUV, HICMA_desc_t *Crk,
                     int rk,
                     int maxrk,
                     double acc);

int HICMA_zgemm_Tile_Async(HICMA_enum transA, HICMA_enum transB,
                           double alpha,
                           HICMA_desc_t *AUV, HICMA_desc_t *Ark,
                           HICMA_desc_t *BUV, HICMA_desc_t *Brk,
                           double beta,
                           HICMA_desc_t *CUV, HICMA_desc_t *Crk,
                           int rk,
                           int maxrk,
                           double acc,
                           HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zgytlr(
        HICMA_enum uplo,
        int M, int N,
        double *AUV,
        double *AD,
        double *Ark,
        int LDA, unsigned long long int seed,
        int maxrank, double tol
);

int HICMA_zgytlr_Tile(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        unsigned long long int seed,
        int maxrank,
        double tol,
        int compress_diag,
        HICMA_desc_t *Dense
);

int HICMA_zgytlr_Tile_Async(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zhagcm_Tile(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank,
        double tol);

int HICMA_zhagcm_Tile_Async(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank, double tol,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zhagdm_Tile(
        HICMA_enum uplo,
        HICMA_desc_t *Dense);

int HICMA_zhagdm_Tile_Async(
        HICMA_enum uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence,
        HICMA_request_t *request);

int HICMA_zhagdmdiag_Tile(
        HICMA_enum uplo,
        HICMA_desc_t *Dense);

int HICMA_zhagdmdiag_Tile_Async(
        HICMA_enum uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence,
        HICMA_request_t *request);

int HICMA_ztrsm_Tile(HICMA_enum side, HICMA_enum uplo,
                     HICMA_enum transA, HICMA_enum diag,
                     double alpha,
                     HICMA_desc_t *AUV,
                     HICMA_desc_t *AD,
                     HICMA_desc_t *Ark,
                     HICMA_desc_t *BUV,
                     HICMA_desc_t *Brk,
                     int rk,
                     int maxrk,
                     double acc);

int HICMA_ztrsm_Tile_Async(HICMA_enum side, HICMA_enum uplo,
                           HICMA_enum transA, HICMA_enum diag,
                           double alpha,
                           HICMA_desc_t *AUV,
                           HICMA_desc_t *AD,
                           HICMA_desc_t *Ark,
                           HICMA_desc_t *BUV,
                           HICMA_desc_t *Brk,
                           int rk,
                           int maxrk,
                           double acc,
                           HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_ztrsmd_Tile(HICMA_enum side, HICMA_enum uplo,
                      HICMA_enum transA, HICMA_enum diag,
                      double alpha,
                      HICMA_desc_t *AUV,
                      HICMA_desc_t *AD,
                      HICMA_desc_t *Ark,
                      HICMA_desc_t *Bdense,
                      int maxrk);

int HICMA_ztrsmd_Tile_Async(HICMA_enum side, HICMA_enum uplo,
                            HICMA_enum transA, HICMA_enum diag,
                            double alpha,
                            HICMA_desc_t *AUV,
                            HICMA_desc_t *AD,
                            HICMA_desc_t *Ark,
                            HICMA_desc_t *Bdense,
                            int maxrk,
                            HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zuncompress(
        HICMA_enum uplo, HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark);

int HICMA_zuncompress_custom_size(HICMA_enum uplo,
                                  HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark,
                                  int numrows_matrix,
                                  int numcolumns_matrix,
                                  int numrows_block,
                                  int numcolumns_block);

int HICMA_zdiag_vec2mat(
        HICMA_desc_t *vec, HICMA_desc_t *mat);

int HICMA_zgenmat_Tile(
        HICMA_desc_t *A);

int HICMA_zgenmat_Tile_Async(
        HICMA_desc_t *A,
        HICMA_sequence_t *sequence,
        HICMA_request_t *request);

int HICMA_zgenrhs_Tile(
        HICMA_desc_t *A);

int HICMA_zgenrhs_Tile_Async(
        HICMA_desc_t *A,
        HICMA_sequence_t *sequence,
        HICMA_request_t *request);

int HICMA_zlacpy(HICMA_enum uplo, int M, int N,
                       HICMA_Complex64_t *A, int LDA,
                       HICMA_Complex64_t *B, int LDB);

int HICMA_zlacpy_Tile(HICMA_enum uplo, HICMA_desc_t *A, HICMA_desc_t *B);

int HICMA_zlacpy_Tile_Async(HICMA_enum uplo, HICMA_desc_t *A, HICMA_desc_t *B,
                                  HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zlaset(HICMA_enum uplo, int M, int N,
                       HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                       HICMA_Complex64_t *A, int LDA);

int HICMA_zlaset_Tile(HICMA_enum uplo,
                            HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                            HICMA_desc_t *A);

int HICMA_zlaset_Tile_Async(HICMA_enum uplo,
                                  HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                                  HICMA_desc_t *A,
                                  HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zgetrf_Tile(HICMA_enum uplo,
                      HICMA_desc_t *AUV,
                      HICMA_desc_t *AD,
                      HICMA_desc_t *Ark,
                      int rk, int maxrk, double acc);

int HICMA_zgetrf_Tile_Async(HICMA_enum uplo,
                            HICMA_desc_t *AUV,
                            HICMA_desc_t *AD,
                            HICMA_desc_t *Ark,
                            int rk, int maxrk, double acc,
                            HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zgenrhs_Tile(HICMA_desc_t *A);

int HICMA_zgenrhs_Tile_Async(HICMA_desc_t *A, HICMA_sequence_t *sequence, HICMA_request_t *request);

int HICMA_zhagdm_Tile(HICMA_enum   uplo,HICMA_desc_t *Dense);

int HICMA_zhagdm_Tile_Async(
        HICMA_enum       uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence,
        HICMA_request_t  *request);

int HICMA_zhagdmdiag_Tile(
        HICMA_enum   uplo,
        HICMA_desc_t *Dense);

int HICMA_zhagdmdiag_Tile_Async(
        HICMA_enum       uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence,
        HICMA_request_t  *request);

int HICMA_zplrnt(int M, int N,
                       HICMA_Complex64_t *A, int LDA,
                       unsigned long long int seed);

int HICMA_zplrnt_Tile(HICMA_desc_t *A,
                            unsigned long long int seed);

int HICMA_zplrnt_Tile_Async(HICMA_desc_t *A,
                                  unsigned long long int seed,
                                  HICMA_sequence_t *sequence,
                                  HICMA_request_t *request);

int HICMA_zLapack_to_Tile( HICMA_Complex64_t *Af77, int LDA, HICMA_desc_t *A );

int HICMA_zTile_to_Lapack( HICMA_desc_t *A, HICMA_Complex64_t *Af77, int LDA );

#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
