/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hicma_z.h
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/
#ifndef _HICMA_Z_H_
#define _HICMA_Z_H_

#undef REAL
#define COMPLEX
#include "morse.h"

#ifdef __cplusplus
extern "C" {
#endif


//not implemented
int HICMA_zpotrf(MORSE_enum uplo, int N, double *A, int LDA);
int HICMA_zpotrf_Tile(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark,
        int rk, int maxrk, double acc
        );
int HICMA_zpotrf_Tile_Async(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark,
        int rk, int maxrk, double acc,
        MORSE_sequence_t *sequence, MORSE_request_t *request );
int HICMA_zgemm_Tile(MORSE_enum transA, MORSE_enum transB,
        double alpha,
        MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        MORSE_desc_t *BUV, MORSE_desc_t *Brk,
        double beta,
        MORSE_desc_t *CUV, MORSE_desc_t *Crk ,
        int rk,
        int maxrk,
        double acc
        );
int HICMA_zgemm_Tile_Async(MORSE_enum transA, MORSE_enum transB,
        double alpha,
        MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        MORSE_desc_t *BUV, MORSE_desc_t *Brk,
        double beta,
        MORSE_desc_t *CUV, MORSE_desc_t *Crk,
        int rk,
        int maxrk,
        double acc ,
        MORSE_sequence_t *sequence, MORSE_request_t *request);
int HICMA_zgytlr(
        MORSE_enum   uplo,
        int M, int N,
        double *AUV,
        double *AD,
        double *Ark,
        int LDA, unsigned long long int seed,
        int maxrank, double tol
        );
int HICMA_zgytlr_Tile(
        MORSE_enum   uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        unsigned long long int seed,
        int maxrank,
        double tol,
        int compress_diag,
        MORSE_desc_t *Dense
        );
int HICMA_zgytlr_Tile_Async(
        MORSE_enum   uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request );
#ifdef __cplusplus
}
#endif

#undef COMPLEX

#endif
