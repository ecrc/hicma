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
#ifndef _COMPUTE_HICMA_Z_H_
#define _COMPUTE_HICMA_Z_H_
/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void hicma_pzpotrf(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark,
        MORSE_sequence_t *sequence, MORSE_request_t *request,
        int rk, int maxrk, double acc);
void hicma_pzgytlr(
        MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request);
void hicma_pzgemm(MORSE_enum transA, MORSE_enum transB,
                  double alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
                  // MORSE_Complex64_t alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
                                           MORSE_desc_t *BUV, MORSE_desc_t *Brk,
                  double beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
                  // MORSE_Complex64_t beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
                  MORSE_sequence_t *sequence, MORSE_request_t *request,
                  int rk, int maxrk, double acc);
#endif
