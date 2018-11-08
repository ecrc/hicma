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
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
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
void hicma_pzhagcm(
        MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank, double tol,
        MORSE_sequence_t *sequence, MORSE_request_t *request );
void hicma_pzhagdm(
        MORSE_enum uplo,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request );
void hicma_pzhagdmdiag(
        MORSE_enum uplo,
        MORSE_desc_t *Dense,
        MORSE_sequence_t *sequence, MORSE_request_t *request );
void hicma_pzgemm(MORSE_enum transA, MORSE_enum transB,
        double alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        // MORSE_Complex64_t alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
        MORSE_desc_t *BUV, MORSE_desc_t *Brk,
        double beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
        // MORSE_Complex64_t beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
        MORSE_sequence_t *sequence, MORSE_request_t *request,
        int rk, int maxrk, double acc);//FIXME put sequence and request at the end
void hicma_pztrsm(MORSE_enum side, MORSE_enum uplo, MORSE_enum trans, MORSE_enum diag,
        double alpha, 
        MORSE_desc_t *AUV, 
        MORSE_desc_t *AD, 
        MORSE_desc_t *Ark, 
        MORSE_desc_t *BUV,
        MORSE_desc_t *Brk,
        int rk,
        int maxrk,
        double acc,
        MORSE_sequence_t *sequence, MORSE_request_t *request);
void hicma_pztrsmd(MORSE_enum side, MORSE_enum uplo, MORSE_enum trans, MORSE_enum diag,
        double alpha, 
        MORSE_desc_t *AUV, 
        MORSE_desc_t *AD, 
        MORSE_desc_t *Ark, 
        MORSE_desc_t *Bdense,
        int maxrk,
        MORSE_sequence_t *sequence, MORSE_request_t *request);
#endif
