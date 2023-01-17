/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 * @file hicma_compute_z.h
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/


#ifndef _COMPUTE_HICMA_Z_H_
#define _COMPUTE_HICMA_Z_H_

#include <control/hicma_descriptor.h>
#include <control/common.h>

/***************************************************************************//**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 **/
void hicma_pzpotrf(HICMA_enum uplo,
                   HICMA_desc_t *AUV, HICMA_desc_t *AD, HICMA_desc_t *Ark,
                   HICMA_sequence_t *sequence, HICMA_request_t *request,
                   int rk, int maxrk, double acc);

void hicma_pzgytlr(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        unsigned long long int seed,
        int maxrank, double tol,
        int compress_diag,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzhagcm(
        HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *Ark,
        int numrows_matrix,
        int numcols_matrix,
        int numrows_block,
        int numcols_block,
        int maxrank, double tol,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzhagdm(
        HICMA_enum uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzhagdmdiag(
        HICMA_enum uplo,
        HICMA_desc_t *Dense,
        HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzgemm(HICMA_enum transA, HICMA_enum transB,
                  double alpha, HICMA_desc_t *AUV, HICMA_desc_t *Ark,
        // HICMA_Complex64_t alpha, HICMA_desc_t *AUV, HICMA_desc_t *Ark,
                  HICMA_desc_t *BUV, HICMA_desc_t *Brk,
                  double beta, HICMA_desc_t *CUV, HICMA_desc_t *Crk,
        // HICMA_Complex64_t beta,  HICMA_desc_t *CUV, HICMA_desc_t *Crk,
                  HICMA_sequence_t *sequence, HICMA_request_t *request,
                  int rk, int maxrk, double acc);//FIXME put sequence and request at the end
void hicma_pztrsm(HICMA_enum side, HICMA_enum uplo, HICMA_enum trans, HICMA_enum diag,
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

void hicma_pztrsmd(HICMA_enum side, HICMA_enum uplo, HICMA_enum trans, HICMA_enum diag,
                   double alpha,
                   HICMA_desc_t *AUV,
                   HICMA_desc_t *AD,
                   HICMA_desc_t *Ark,
                   HICMA_desc_t *Bdense,
                   int maxrk,
                   HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzlacpy(HICMA_enum uplo, HICMA_desc_t *A, HICMA_desc_t *B,
                         HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzlaset(HICMA_enum uplo,
                         HICMA_Complex64_t alpha, HICMA_Complex64_t beta,
                         HICMA_desc_t *A,
                         HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzplrnt(HICMA_desc_t *A, unsigned long long int seed,
                         HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzgenrhs(HICMA_desc_t *A, HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzgenmat(HICMA_desc_t *A, HICMA_sequence_t *sequence, HICMA_request_t *request);

void hicma_pzgetrf(HICMA_enum uplo,
                   HICMA_desc_t *AUV,
                   HICMA_desc_t *AD,
                   HICMA_desc_t *Ark,
                   HICMA_sequence_t *sequence, HICMA_request_t *request,
                   int rk, int maxrk, double acc);

/**
 *
 * @file compute_z.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon computational functions header
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for HICMA 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */
/**
 *  LAPACK/Tile Descriptor accesses
 */
#define HicmaDescInput  1
#define HicmaDescOutput 2
#define HicmaDescInout  (HicmaDescInput | HicmaDescOutput)

/**
 *  Macro for matrix conversion / Lapack interface
 */
#define hicma_zdesc_alloc_diag( descA, mb, nb, lm, ln, i, j, m, n, p, q) \
    descA = hicma_desc_init_diag(                                       \
        HicmaComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), p, q);                            \
    hicma_desc_mat_alloc( &(descA) );                                   \
    HICMA_RUNTIME_desc_create( &(descA) );

#define hicma_zdesc_alloc( descA, mb, nb, lm, ln, i, j, m, n, free)     \
    descA = hicma_desc_init(                                            \
        HicmaComplexDouble, (mb), (nb), ((mb)*(nb)),                    \
        (m), (n), (i), (j), (m), (n), 1, 1);                            \
    if ( hicma_desc_mat_alloc( &(descA) ) ) {                           \
        hicma_error( __func__, "hicma_desc_mat_alloc() failed");        \
        {free;};                                                        \
        return HICMA_ERR_OUT_OF_RESOURCES;                              \
    }                                                                   \
    HICMA_RUNTIME_desc_create( &(descA) );

/**
 *  Declarations of parallel functions (dynamic scheduling) - alphabetical order
 */
void hicma_pzlacpy(HICMA_enum uplo, HICMA_desc_t *A, HICMA_desc_t *B, HICMA_sequence_t *sequence, HICMA_request_t *request);

/**
 * @brief Internal function to convert the lapack format to tile format in
 * LAPACK interface calls
 */
static inline int
hicma_zlap2tile( HICMA_context_t *hicma,
                 HICMA_desc_t *descAl, HICMA_desc_t *descAt,
                 HICMA_enum mode, HICMA_enum uplo,
                 HICMA_Complex64_t *A, int mb, int nb, int lm, int ln, int m, int n,
                 HICMA_sequence_t *seq, HICMA_request_t *req )
{
    /* Initialize the Lapack descriptor */
    *descAl = hicma_desc_init_user( HicmaComplexDouble, mb, nb, (mb)*(nb),
                                    lm, ln, 0, 0, m, n, 1, 1,
                                    hicma_getaddr_cm, hicma_getblkldd_cm, NULL  );
    descAl->mat = A;
    descAl->styp = HicmaCM;

    /* Initialize the tile descriptor */
    *descAt = hicma_desc_init( HicmaComplexDouble, mb, nb, (mb)*(nb),
                               lm, ln, 0, 0, m, n, 1, 1 );

    if ( HICMA_TRANSLATION == HICMA_OUTOFPLACE ) {
        if ( hicma_desc_mat_alloc( descAt ) ) {
            hicma_error( "hicma_zlap2tile", "hicma_desc_mat_alloc() failed");
            return HICMA_ERR_OUT_OF_RESOURCES;
        }

        HICMA_RUNTIME_desc_create( descAl );
        HICMA_RUNTIME_desc_create( descAt );

        if ( mode & HicmaDescInput ) {
            hicma_pzlacpy( uplo, descAl, descAt, seq, req );
        }

    }
    else {
        hicma_fatal_error( "hicma_zlap2tile", "INPLACE translation not supported yet");
        descAt->mat = A;

        HICMA_RUNTIME_desc_create( descAl );
        HICMA_RUNTIME_desc_create( descAt );

        if ( mode & HicmaDescInput ) {
            /* HICMA_zgecfi_Async( lm, ln, A, HicmaCM, mb, nb, */
            /*                     HicmaCCRB, mb, nb, seq, req ); */
        }
        return HICMA_ERR_NOT_SUPPORTED;
    }
    return HICMA_SUCCESS;
}

/**
 * @brief Internal function to convert back the tile format to the lapack format
 * in LAPACK interface calls
 */
static inline int
hicma_ztile2lap( HICMA_context_t *hicma, HICMA_desc_t *descAl, HICMA_desc_t *descAt,
                 HICMA_enum mode, HICMA_enum uplo, HICMA_sequence_t *seq, HICMA_request_t *req )
{
    if ( HICMA_TRANSLATION == HICMA_OUTOFPLACE ) {
        if ( mode & HicmaDescOutput ) {
            hicma_pzlacpy( uplo, descAt, descAl, seq, req );
        }
    }
    else {
        hicma_fatal_error( "hicma_ztile2lap", "INPLACE translation not supported yet");
        if ( mode & HicmaDescOutput ) {
            /* HICMA_zgecfi_Async( descAl->lm, descAl->ln, descAl->mat, */
            /*                     HicmaCCRB, descAl->mb, descAl->nb,   */
            /*                     HicmaCM, descAl->mb, descAl->nb, seq, req ); */
        }
        return HICMA_ERR_NOT_SUPPORTED;
    }
    HICMA_RUNTIME_desc_flush( descAl, seq );
    HICMA_RUNTIME_desc_flush( descAt, seq );

    return HICMA_SUCCESS;
}

/**
 * @brief Internal function to cleanup the temporary data from the layout
 * conversions in LAPACK interface calls
 */
static inline void
hicma_ztile2lap_cleanup( HICMA_context_t *hicma, HICMA_desc_t *descAl, HICMA_desc_t *descAt )
{
    if ( HICMA_TRANSLATION == HICMA_OUTOFPLACE ) {
        hicma_desc_mat_free( descAt );
    }
    HICMA_RUNTIME_desc_destroy( descAl );
    HICMA_RUNTIME_desc_destroy( descAt );
}

#endif
