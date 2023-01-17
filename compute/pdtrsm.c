/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pdtrsm.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

/**
 *
 * file pdtrsm.c
 *
 * MORSE auxiliary routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * @author Jakub Kurzak
 * @author Hatem Ltaief
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 **/

#include <control/common.h>
#include "hicma.h"
#include <hicma_runtime_d.h>
#include <stdio.h>
#include <assert.h>
#include <cblas.h>

#define A(m, n) A,  m,  n
#define B(m, n) B,  m,  n
int pdtrsm_enable_dense = 0;

/***************************************************************************//**
 *  Parallel tile triangular solve - dynamic scheduling
 **/
void hicma_pdtrsm(HICMA_enum side, HICMA_enum uplo, HICMA_enum trans, HICMA_enum diag,
                  double alpha,
                  HICMA_desc_t *AUV,
                  HICMA_desc_t *AD,
                  HICMA_desc_t *Ark,
                  HICMA_desc_t *BUV,
                  HICMA_desc_t *Brk,
                  int rk,
                  int maxrk,
                  double acc,
                  HICMA_sequence_t *sequence, HICMA_request_t *request) {
    if (HICMA_get_print_index() == 1) {
        printf("%d:%s rk:%d maxrk:%d acc:%e alpha:%e\n",
               __LINE__, __func__,
               rk, maxrk, acc, alpha);
    }
    HICMA_desc_t *A = AUV;
    HICMA_desc_t *B = BUV;

    HICMA_context_t *hicma;
    HICMA_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    double done = (double) 1.0;
    double mdone = (double) -1.0;
    double minvalpha = (double) -1.0 / alpha;
    double lalpha;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);
    size_t ws_host = 0;
    size_t ws_worker = 0;
    ws_worker =  //FIXME tentative size. Find exact size. I think syrk uses less memory
            2 * AD->mb * 2 *
            maxrk   // for copying CU and CV into temporary buffer instead of using CUV itself. There is 2*maxrk because these buffers will be used to put two U's side by side
            + 2 * AD->mb       // qrtauA qrtauB
            + maxrk * maxrk    // qrb_aubut  AcolBcolT
            + 2 * AD->mb * maxrk // newU newV
            + (2 * maxrk) * (2 * maxrk)    // svd_rA     _rA
            + (2 * maxrk)          // sigma
            ;
    ws_worker *= sizeof(double); //FIXME use HICMA_Complex64_t
    //FIXME add ws_worker and ws_host calculation from compute/pzgeqrf.c when GPU/MAGMA is supported
    HICMA_RUNTIME_options_ws_alloc(&options, ws_worker, ws_host);
    /*
     *  HicmaLeft / HicmaUpper / HicmaNoTrans
     */
    if (side == HicmaLeft) {
        if (uplo == HicmaUpper) {
            assert("Not implemented yet" == 0);
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m - (B->mt - 1) * B->mb : B->mb;
                    ldak = BLKLDD(A, B->mt - 1 - k);
                    ldbk = BLKLDD(B, B->mt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempkm, tempnn, A->mb,
//                                lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,  /* lda * tempkm */
//                                B(B->mt - 1 - k, n), ldbk); /* ldb * tempnn */
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,
                                               B(B->mt - 1 - k, n), ldbk);

                    }
                    //HICMA_RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-k) );
                    HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, B->mt - 1 - k));

                    for (m = k + 1; m < B->mt; m++) {
                        ldam = BLKLDD(A, B->mt - 1 - m);
                        ldbm = BLKLDD(B, B->mt - 1 - m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    B->mb, tempnn, tempkm, A->mb,
//                                    mdone, A(B->mt - 1 - m, B->mt - 1 - k), ldam,
//                                    B(B->mt - 1 - k, n), ldbk,
//                                    lalpha, B(B->mt - 1 - m, n), ldbm);


                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             B->mb, tempnn, tempkm, mdone,
                                             A(B->mt - 1 - m, B->mt - 1 - k), ldam,
                                             B(B->mt - 1 - k, n), ldbk,
                                             lalpha, B(B->mt - 1 - m, n), ldbm);

                        }
                        //HICMA_RUNTIME_data_flush( sequence, A(B->mt-1-m, B->mt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - m, B->mt - 1 - k));
                    }
                    for (n = 0; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, B(B->mt-1-k, n) );
                        HICMA_RUNTIME_data_flush(sequence, B(B->mt - 1 - k, n));
                    }
                }
            }
                /*
                 *  HicmaLeft / HicmaUpper / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt - 1 ? B->m - k * B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempkm, tempnn, A->mb,
//                                lalpha, A(k, k), ldak,
//                                B(k, n), ldbk);
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, A(k, k), ldak,
                                               B(k, n), ldbk);

                    }
                    //HICMA_RUNTIME_data_flush( sequence, A(k, k) );
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = k + 1; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    trans, HicmaNoTrans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    mdone, A(k, m), ldak,
//                                    B(k, n), ldbk,
//                                    lalpha, B(m, n), ldbm);

                            HICMA_TASK_dgemm(&options, trans, HicmaNoTrans,
                                             tempmm, tempnn, B->mb, mdone,
                                             A(k, m), ldak,
                                             B(k, n), ldbk,
                                             lalpha, B(m, n), ldbm);

                        }
                        //HICMA_RUNTIME_data_flush( sequence, A(k, m) );
                        HICMA_RUNTIME_data_flush(sequence, A(k, m));
                    }
                    for (n = 0; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, B(k, n) );
                        HICMA_RUNTIME_data_flush(sequence, B(k, n));
                    }

                }
            }
        }
            /*
             *  HicmaLeft / HicmaLower / HicmaNoTrans
             */
        else {
            if (trans == HicmaNoTrans) {
                //@1
                //printf("%s %d Left Lower Notrans\n", __FILE__, __LINE__);
                for (k = 0; k < B->mt; k++) {
                    int ldbkuv = BLKLDD(BUV, k);
                    int ldakd = BLKLDD(AD, k);
                    tempkm = k == B->mt - 1 ? B->m - k * B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                        if (pdtrsm_enable_dense) {

//                            HICMA_TASK_hcore_dtrsm(
//                                    &options,
//                                    side, uplo, trans, diag,
//                                    tempkm, tempnn, A->mb,
//                                    lalpha, A(k, k), ldak,
//                                    B(k, n), ldbk);

                            HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                                   tempkm, tempnn, lalpha, A(k, k), ldak,
                                                   B(k, n), ldbk);
                        } else {
                            HICMA_TASK_hcore_dtrsm(
                                    &options,
                                    //HicmaLeft, HicmaLower, HicmaNoTrans, HicmaNonUnit,
                                    side, uplo, trans, diag,
                                    tempkm, //FIXME must be number of rows of the diagonal block
                                    lalpha, AD, k,
                                    0, // I assume that only diags are stored
                                    ldakd,
                                    BUV, k, n,
                                    ldbkuv,
                                    Brk);
                        }
                    }
                    //HICMA_RUNTIME_data_flush( sequence, A(k, k) );
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = k + 1; m < B->mt; m++) {
                        int ldamuv = BLKLDD(AUV, m);
                        int ldbmuv = BLKLDD(BUV, m);
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldam = BLKLDD(A, m);
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            if (pdtrsm_enable_dense) {

//                                HICMA_TASK_dgemm(
//                                        &options,
//                                        HicmaNoTrans, HicmaNoTrans,
//                                        tempmm, tempnn, B->mb, A->mb,
//                                        mdone, A(m, k), ldam,
//                                        B(k, n), ldbk,
//                                        lalpha, B(m, n), ldbm);

                                HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                                 tempmm, tempnn, B->mb, mdone, A(m, k), ldam, B(k, n), ldbk,
                                                 lalpha, B(m, n), ldbm);

                            } else {
                                HICMA_TASK_hcore_dgemm(
                                        &options,
                                        HicmaNoTrans, HicmaNoTrans,
                                        tempmm, //TODO tempmmuv,
                                        tempnn, //TODO tempmmuv,
                                        mdone,
                                        AUV, Ark, m, k, ldamuv,
                                        BUV, Brk, k, n, ldbkuv,
                                        lalpha,
                                        BUV, Brk, m, n, ldbmuv,
                                        rk, maxrk, acc);
                            }
                        }
                        //HICMA_RUNTIME_data_flush( sequence, A(m, k) );
                        HICMA_RUNTIME_data_flush(sequence, A(m, k));
                    }
                    for (n = 0; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, B(k, n) );
                        HICMA_RUNTIME_data_flush(sequence, B(k, n));
                    }
                }
            }
                /*
                 *  HicmaLeft / HicmaLower / Hicma[Conj]Trans
                 */
            else {
                //@2
                //printf("%s %d Left Lower Trans\n", __FILE__, __LINE__);
                for (k = 0; k < B->mt; k++) {
                    int ldakuv = BLKLDD(AUV, B->mt - 1 - k);
                    int ldbkuv = BLKLDD(BUV, B->mt - 1 - k);
                    int ldakd = BLKLDD(AD, B->mt - 1 - k);
                    tempkm = k == 0 ? B->m - (B->mt - 1) * B->mb : B->mb;
                    ldak = BLKLDD(A, B->mt - 1 - k);
                    ldbk = BLKLDD(B, B->mt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                        if (pdtrsm_enable_dense) {
//                            HICMA_TASK_hcore_dtrsm(
//                                    &options,
//                                    side, uplo, trans, diag,
//                                    tempkm, tempnn, A->mb,
//                                    lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,
//                                    B(B->mt - 1 - k, n), ldbk);
                            HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                                   tempkm, tempnn, lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,
                                                   B(B->mt - 1 - k, n), ldbk);

                        } else {
                            HICMA_TASK_hcore_dtrsm(
                                    &options,
                                    //HicmaLeft, HicmaLower, HicmaNoTrans, HicmaNonUnit,
                                    side, uplo, trans, diag,
                                    tempkm, //FIXME must be number of rows of the diagonal block
                                    lalpha, AD, B->mt - 1 - k,
                                    0, // I assume that only diags are stored
                                    ldakd,
                                    BUV, k, n,
                                    ldbkuv,
                                    Brk);
                        }
                    }
                    //HICMA_RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-k) );
                    HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, B->mt - 1 - k));
                    for (m = k + 1; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, B->mt - 1 - m);
                        int ldbmuv = BLKLDD(BUV, B->mt - 1 - m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            if (pdtrsm_enable_dense) {
//                                HICMA_TASK_dgemm(
//                                        &options,
//                                        trans, HicmaNoTrans,
//                                        B->mb, tempnn, tempkm, A->mb,
//                                        mdone, A(B->mt - 1 - k, B->mt - 1 - m), ldak,
//                                        B(B->mt - 1 - k, n), ldbk,
//                                        lalpha, B(B->mt - 1 - m, n), ldbm);

                                HICMA_TASK_dgemm(&options, trans, HicmaNoTrans,
                                                 B->mb, tempnn, tempkm, mdone,
                                                 A(B->mt - 1 - k, B->mt - 1 - m), ldak,
                                                 B(B->mt - 1 - k, n), ldbk,
                                                 lalpha, B(B->mt - 1 - m, n), ldbm);

                            } else {
                                HICMA_TASK_hcore_dgemm(
                                        &options,
                                        trans, HicmaNoTrans,
                                        B->mb, //TODO tempmmuv,
                                        tempnn, //TODO tempmmuv,
                                        mdone,
                                        AUV, Ark, B->mt - 1 - k, B->mt - 1 - m, ldakuv,
                                        BUV, Brk, B->mt - 1 - k, n, ldbkuv,
                                        lalpha,
                                        BUV, Brk, B->mt - 1 - m, n, ldbmuv,
                                        rk, maxrk, acc);
                            }
                        }
                        //HICMA_RUNTIME_data_flush( sequence, A(B->mt-1-k, B->mt-1-m) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, B->mt - 1 - m));
                    }
                    for (n = 0; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, B(B->mt-1-k, n) );
                        HICMA_RUNTIME_data_flush(sequence, B(B->mt - 1 - k, n));
                    }
                }
            }
        }
    }
        /*
         *  HicmaRight / HicmaUpper / HicmaNoTrans
         */
    else {
        assert("Not implemented yet" == 0);
        if (uplo == HicmaUpper) {
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt - 1 ? B->n - k * B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : done;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                lalpha, A(k, k), ldak,  /* lda * tempkn */
//                                B(m, k), ldbm); /* ldb * tempkn */
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, A(k, k), ldak,
                                               B(m, k), ldbm);

                    }
                    //HICMA_RUNTIME_data_flush( sequence, A(k, k) );
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = k + 1; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    mdone, B(m, k), ldbm,  /* ldb * B->mb   */
//                                    A(k, n), ldak,  /* lda * tempnn */
//                                    lalpha, B(m, n), ldbm); /* ldb * tempnn */


                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             tempmm, tempnn, B->mb, mdone,
                                             B(m, k), ldbm,
                                             A(k, n), ldak,
                                             lalpha, B(m, n), ldbm);

                        }
                        //HICMA_RUNTIME_data_flush( sequence, B(m, k) );
                        HICMA_RUNTIME_data_flush(sequence, B(m, k));
                    }
                    for (n = k + 1; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, A(k, n) );
                        HICMA_RUNTIME_data_flush(sequence, A(k, n));
                    }
//                    free(A_array);
//                    free(B_array);
//                    free(C_array);
                }
            }
                /*
                 *  HicmaRight / HicmaUpper / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n - (B->nt - 1) * B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt - 1 - k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                alpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,  /* lda * tempkn */
//                                B(m, B->nt - 1 - k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, alpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,
                                               B(m, B->nt - 1 - k), ldbm);

                        //HICMA_RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - k));

                        for (n = k + 1; n < B->nt; n++) {
                            ldan = BLKLDD(A, B->nt - 1 - n);
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, trans,
//                                    tempmm, B->nb, tempkn, A->mb,
//                                    minvalpha, B(m, B->nt - 1 - k), ldbm,  /* ldb  * tempkn */
//                                    A(B->nt - 1 - n, B->nt - 1 - k), ldan, /* A->mb * tempkn (Never last row) */
//                                    done, B(m, B->nt - 1 - n), ldbm); /* ldb  * B->nb   */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, trans,
                                             tempmm, B->nb, tempkn, minvalpha,
                                             B(m, B->nt - 1 - k), ldbm,
                                             A(B->nt - 1 - n, B->nt - 1 - k), ldan,
                                             done, B(m, B->nt - 1 - n), ldbm);

                        }
                        //HICMA_RUNTIME_data_flush( sequence, B(m,        B->nt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, B(m, B->nt - 1 - k));
                    }
                    for (n = k + 1; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, A(B->nt-1-n, B->nt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - n, B->nt - 1 - k));
                    }
                }
            }
        }
            /*
             *  HicmaRight / HicmaLower / HicmaNoTrans
             */
        else {
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n - (B->nt - 1) * B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                lalpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,  /* lda * tempkn */
//                                B(m, B->nt - 1 - k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, lalpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,
                                               B(m, B->nt - 1 - k), ldbm);

                        //HICMA_RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - k));

                        for (n = k + 1; n < B->nt; n++) {
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    tempmm, B->nb, tempkn, A->mb,
//                                    mdone, B(m, B->nt - 1 - k), ldbm,  /* ldb * tempkn */
//                                    A(B->nt - 1 - k, B->nt - 1 - n), ldak,  /* lda * B->nb   */
//                                    lalpha, B(m, B->nt - 1 - n), ldbm); /* ldb * B->nb   */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             tempmm, B->nb, tempkn, mdone,
                                             B(m, B->nt - 1 - k), ldbm,
                                             A(B->nt - 1 - k, B->nt - 1 - n), ldak,
                                             lalpha, B(m, B->nt - 1 - n), ldbm);

                        }

                        //HICMA_RUNTIME_data_flush( sequence, B(m,        B->nt-1-k) );
                        HICMA_RUNTIME_data_flush(sequence, B(m, B->nt - 1 - k));

                    }
                    for (n = k + 1; n < B->nt; n++) {
                        //HICMA_RUNTIME_data_flush( sequence, A(B->nt-1-k, B->nt-1-n) );
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - n));
                    }
                }
            }
                /*
                 *  HicmaRight / HicmaLower / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt - 1 ? B->n - k * B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                alpha, A(k, k), ldak,  /* lda * tempkn */
//                                B(m, k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, alpha, A(k, k), ldak,
                                               B(m, k), ldbm);

                        //HICMA_RUNTIME_data_flush( sequence, A(k, k) );
                        HICMA_RUNTIME_data_flush(sequence, A(k, k));

                        for (n = k + 1; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            ldan = BLKLDD(A, n);
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, trans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    minvalpha, B(m, k), ldbm,  /* ldb  * tempkn */
//                                    A(n, k), ldan, /* ldan * tempkn */
//                                    done, B(m, n), ldbm); /* ldb  * tempnn */
//
                            HICMA_TASK_dgemm(&options, HicmaNoTrans, trans,
                                             tempmm, tempnn, B->mb, minvalpha,
                                             B(m, k), ldbm,
                                             A(n, k), ldan,
                                             done, B(m, n), ldbm);
                        }
                        HICMA_RUNTIME_data_flush(sequence, B(m, k));

                    }
                    for (n = k + 1; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, A(n, k));
                    }

                }
            }
        }
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicma);
}

void hicma_pdtrsmd(HICMA_enum side, HICMA_enum uplo, HICMA_enum trans, HICMA_enum diag,
                   double alpha,
                   HICMA_desc_t *AUV,
                   HICMA_desc_t *AD,
                   HICMA_desc_t *Ark,
                   HICMA_desc_t *Bdense,
                   int maxrk,
                   HICMA_sequence_t *sequence, HICMA_request_t *request) {
    if (HICMA_get_print_index() == 1) {
        printf("%d:%s maxrk:%d alpha:%e\n",
               __LINE__, __func__,
               maxrk, alpha);
    }
    HICMA_desc_t *A = AUV;
    HICMA_desc_t *B = Bdense;

    HICMA_context_t *hicma;
    HICMA_option_t options;

    int k, m, n;
    int ldak, ldam, ldan, ldbk, ldbm;
    int tempkm, tempkn, tempmm, tempnn;

    double done = (double) 1.0;
    double mdone = (double) -1.0;
    double minvalpha = (double) -1.0 / alpha;
    double lalpha;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);
    size_t ws_host = 0;
    size_t ws_worker = 0;
    ws_worker =
            +AD->mb * maxrk // temporary space for performing AV*B in CD+=AU*(AV*B)
            ;
    ws_worker *= sizeof(double); //FIXME use HICMA_Complex64_t
    //FIXME add ws_worker and ws_host calculation from compute/pzgeqrf.c when GPU/MAGMA is supported
    HICMA_RUNTIME_options_ws_alloc(&options, ws_worker, ws_host);
    /*
     *  HicmaLeft / HicmaUpper / HicmaNoTrans
     */
    if (side == HicmaLeft) {
        if (uplo == HicmaUpper) {
            assert("Not implemented yet" == 0);
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == 0 ? B->m - (B->mt - 1) * B->mb : B->mb;
                    ldak = BLKLDD(A, B->mt - 1 - k);
                    ldbk = BLKLDD(B, B->mt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempkm, tempnn, A->mb,
//                                lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,  /* lda * tempkm */
//                                B(B->mt - 1 - k, n), ldbk); /* ldb * tempnn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, A(B->mt - 1 - k, B->mt - 1 - k), ldak,
                                               B(B->mt - 1 - k, n), ldbk);

                    }
                    HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, B->mt - 1 - k));

                    for (m = k + 1; m < B->mt; m++) {
                        ldam = BLKLDD(A, B->mt - 1 - m);
                        ldbm = BLKLDD(B, B->mt - 1 - m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    B->mb, tempnn, tempkm, A->mb,
//                                    mdone, A(B->mt - 1 - m, B->mt - 1 - k), ldam,
//                                    B(B->mt - 1 - k, n), ldbk,
//                                    lalpha, B(B->mt - 1 - m, n), ldbm);

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             B->mb, tempnn, tempkm, mdone,
                                             A(B->mt - 1 - m, B->mt - 1 - k), ldam,
                                             B(B->mt - 1 - k, n), ldbk,
                                             lalpha, B(B->mt - 1 - m, n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - m, B->mt - 1 - k));
                    }
                    for (n = 0; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, B(B->mt - 1 - k, n));
                    }
                }
            }
                /*
                 *  HicmaLeft / HicmaUpper / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->mt; k++) {
                    tempkm = k == B->mt - 1 ? B->m - k * B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempkm, tempnn, A->mb,
//                                lalpha, A(k, k), ldak,
//                                B(k, n), ldbk);
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, A(k, k), ldak,
                                               B(k, n), ldbk);

                    }
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = k + 1; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    trans, HicmaNoTrans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    mdone, A(k, m), ldak,
//                                    B(k, n), ldbk,
//                                    lalpha, B(m, n), ldbm);

                            HICMA_TASK_dgemm(&options, trans, HicmaNoTrans,
                                             tempmm, tempnn, B->mb, mdone,
                                             A(k, m), ldak,
                                             B(k, n), ldbk,
                                             lalpha, B(m, n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, A(k, m));
                    }
                    for (n = 0; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, B(k, n));
                    }
                }
            }
        }
            /*
             *  HicmaLeft / HicmaLower / HicmaNoTrans
             */
        else {
            if (trans == HicmaNoTrans) {
                //@1
                //printf("%s %d Left Lower Notrans\n", __FILE__, __LINE__);
                for (k = 0; k < B->mt; k++) {
                    int ldbkd = BLKLDD(Bdense, k);
                    int ldakd = BLKLDD(AD, k);
                    tempkm = k == B->mt - 1 ? B->m - k * B->mb : B->mb;
                    ldak = BLKLDD(A, k);
                    ldbk = BLKLDD(B, k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempkm,
//                                tempnn,
//                                A->mb,
//                                lalpha, AD, k, 0, ldakd,
//                                Bdense, k, n, ldbkd);
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, AD, k, 0, ldakd,
                                               Bdense, k, n, ldbkd);

                    }
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = k + 1; m < B->mt; m++) {
                        int ldamuv = BLKLDD(AUV, m);
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldam = BLKLDD(A, m);
                        int ldbmd = BLKLDD(Bdense, m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            if (pdtrsm_enable_dense) {

//                                HICMA_TASK_dgemm(
//                                        &options,
//                                        HicmaNoTrans, HicmaNoTrans,
//                                        tempmm, tempnn, B->mb, A->mb,
//                                        mdone, A(m, k), ldam,
//                                        B(k, n), ldbk,
//                                        lalpha, B(m, n), ldbm);

                                HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                                 tempmm, tempnn, B->mb, mdone,
                                                 A(m, k), ldam,
                                                 B(k, n), ldbk,
                                                 lalpha, B(m, n), ldbm);

                            } else {
                                /*printf("(%d,%d,%d): (%d,%d [%d])=(%d,%d [%d])x(%d,%d [%d]) trans:%d tempmm:%d tempnn:%d alpha:%g\n", k, m, n, */
                                /*m, n, ldbmd,*/
                                /*m, k, ldamuv,*/
                                /*k, n, ldbkd,*/
                                /*trans,*/
                                /*tempmm, tempnn, lalpha*/
                                /*);*/
                                HICMA_TASK_dgemm_bdcd(
                                        &options,
                                        HicmaNoTrans, HicmaNoTrans,
                                        tempmm,
                                        tempnn,
                                        mdone,
                                        AUV, Ark, m, k, ldamuv,
                                        Bdense, k, n, ldbkd,
                                        lalpha,
                                        Bdense, m, n, ldbmd);
                            }
                        }
                        HICMA_RUNTIME_data_flush(sequence, A(m, 0));
                    }
                    for (n = 0; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, B(k, n));
                    }
                }
            }
                /*
                 *  HicmaLeft / HicmaLower / Hicma[Conj]Trans
                 */
            else {
                //@2
                if (0)
                    printf("%s %d Left Lower Trans B->m,n:%d,%d B->mt,nt:%d,%d\n", __FILE__, __LINE__, B->m, B->n,
                           B->mt, B->nt);
                for (k = 0; k < B->mt; k++) {
                    int ldakuv = BLKLDD(AUV, B->mt - 1 - k);
                    int ldakd = BLKLDD(AD, B->mt - 1 - k);
                    int ldbkd = BLKLDD(Bdense, B->mt - 1 - k);
                    tempkm = k == 0 ? B->m - (B->mt - 1) * B->mb : B->mb;
                    ldak = BLKLDD(AD, B->mt - 1 - k);
                    ldbk = BLKLDD(B, B->mt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (n = 0; n < B->nt; n++) {
                        tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                        /*printf("chamtrsm: (%d,%d) A(%d,%d [%d]) B(%d,%d [%d])\n",*/
                        /*k, n,*/
                        /*B->mt-1-k, 0, ldak,*/
                        /*B->mt-1-k, n, ldbk*/
                        /*);*/

//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo,
//                                trans,
//                                diag,
//                                tempkm, tempnn, A->mb,
//                                lalpha, AD, B->mt - 1 - k, 0, ldak,
//                                Bdense, B->mt - 1 - k, n, ldbk);
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempkm, tempnn, lalpha, AD, B->mt - 1 - k, 0, ldak,
                                               Bdense, B->mt - 1 - k, n, ldbk);

                    }
                    HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, B->mt - 1 - k));

                    for (m = k + 1; m < B->mt; m++) {
                        tempmm = 0;
                        if (B->mt - 1 - k == B->mt - 1) {
                            if (AD->m % AUV->mb == 0) {
                                tempmm = AUV->mb;
                            } else {
                                tempmm = AD->m % AUV->mb;
                            }
                        } else {
                            tempmm = AUV->mb;
                        }
                        ldbm = BLKLDD(B, B->mt - 1 - m);
                        for (n = 0; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            if (pdtrsm_enable_dense) {
//                                HICMA_TASK_dgemm(
//                                        &options,
//                                        trans, HicmaNoTrans,
//                                        B->mb, tempnn, tempkm, A->mb,
//                                        mdone, A(B->mt - 1 - k, B->mt - 1 - m), ldak,
//                                        B(B->mt - 1 - k, n), ldbk,
//                                        lalpha, B(B->mt - 1 - m, n), ldbm);

                                HICMA_TASK_dgemm(&options, trans, HicmaNoTrans,
                                                 B->mb, tempnn, tempkm, mdone,
                                                 A(B->mt - 1 - k, B->mt - 1 - m), ldak,
                                                 B(B->mt - 1 - k, n), ldbk,
                                                 lalpha, B(B->mt - 1 - m, n), ldbm);

                            } else {
                                /*printf("(%d,%d,%d): (%d,%d [%d])=(%d,%d [%d])x(%d,%d [%d]) trans:%d tempmm:%d tempnn:%d alpha:%g\n", k, m, n, */
                                /*B->mt-1-m, n, ldbm,*/
                                /*B->mt-1-k, B->mt-1-m, ldakuv,*/
                                /*B->mt-1-k, n, ldbk,*/
                                /*trans,*/
                                /*tempmm, tempnn, lalpha);*/
                                HICMA_TASK_dgemm_bdcd(
                                        &options,
                                        trans, HicmaNoTrans,
                                        tempmm,
                                        tempnn,
                                        mdone,
                                        AUV, Ark, B->mt - 1 - k, B->mt - 1 - m, ldakuv,
                                        Bdense, B->mt - 1 - k, n, ldbk,
                                        lalpha,
                                        Bdense, B->mt - 1 - m, n, ldbm);
                            }
                        }
                        HICMA_RUNTIME_data_flush(sequence, A(B->mt - 1 - k, 0));
                    }
                    for (n = 0; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, B(B->mt - 1 - k, n));
                    }
                }
            }
        }
    }
        /*
         *  HicmaRight / HicmaUpper / HicmaNoTrans
         */
    else {
        assert("Not implemented yet" == 0);
        if (uplo == HicmaUpper) {
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt - 1 ? B->n - k * B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    lalpha = k == 0 ? alpha : done;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                lalpha, A(k, k), ldak,  /* lda * tempkn */
//                                B(m, k), ldbm); /* ldb * tempkn */
                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, lalpha, A(k, k), ldak,
                                               B(m, k), ldbm);

                    }
                    HICMA_RUNTIME_data_flush(sequence, A(k, k));

                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
                        for (n = k + 1; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    mdone, B(m, k), ldbm,  /* ldb * B->mb   */
//                                    A(k, n), ldak,  /* lda * tempnn */
//                                    lalpha, B(m, n), ldbm); /* ldb * tempnn */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             tempmm, tempnn, B->mb, mdone,
                                             B(m, k), ldbm,
                                             A(k, n), ldak,
                                             lalpha, B(m, n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, B(m, k));
                    }

                    for (n = k + 1; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, A(k, n));
                    }
                }
            }
                /*
                 *  HicmaRight / HicmaUpper / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n - (B->nt - 1) * B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt - 1 - k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                alpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,  /* lda * tempkn */
//                                B(m, B->nt - 1 - k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, alpha, A(B->nt - 1 - k, B->nt - 1 - k),
                                               ldak, B(m, B->nt - 1 - k), ldbm);

                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - k));

                        for (n = k + 1; n < B->nt; n++) {
                            ldan = BLKLDD(A, B->nt - 1 - n);
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, trans,
//                                    tempmm, B->nb, tempkn, A->mb,
//                                    minvalpha, B(m, B->nt - 1 - k), ldbm,  /* ldb  * tempkn */
//                                    A(B->nt - 1 - n, B->nt - 1 - k), ldan, /* A->mb * tempkn (Never last row) */
//                                    done, B(m, B->nt - 1 - n), ldbm); /* ldb  * B->nb   */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, trans,
                                             tempmm, B->nb, tempkn, minvalpha,
                                             B(m, B->nt - 1 - k), ldbm,
                                             A(B->nt - 1 - n, B->nt - 1 - k), ldan,
                                             done, B(m, B->nt - 1 - n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, B(m, B->nt - 1 - k));
                    }
                    for (n = k + 1; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - n, B->nt - 1 - k));
                    }
                }
            }
        }
            /*
             *  HicmaRight / HicmaLower / HicmaNoTrans
             */
        else {
            if (trans == HicmaNoTrans) {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == 0 ? B->n - (B->nt - 1) * B->nb : B->nb;
                    ldak = BLKLDD(A, B->nt - 1 - k);
                    lalpha = k == 0 ? alpha : done;
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                lalpha, A(B->nt - 1 - k, B->nt - 1 - k), ldak,  /* lda * tempkn */
//                                B(m, B->nt - 1 - k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, lalpha, A(B->nt - 1 - k, B->nt - 1 - k),
                                               ldak, B(m, B->nt - 1 - k), ldbm);


                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - k));

                        for (n = k + 1; n < B->nt; n++) {
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, HicmaNoTrans,
//                                    tempmm, B->nb, tempkn, A->mb,
//                                    mdone, B(m, B->nt - 1 - k), ldbm,  /* ldb * tempkn */
//                                    A(B->nt - 1 - k, B->nt - 1 - n), ldak,  /* lda * B->nb   */
//                                    lalpha, B(m, B->nt - 1 - n), ldbm); /* ldb * B->nb   */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, HicmaNoTrans,
                                             tempmm, B->nb, tempkn, mdone,
                                             B(m, B->nt - 1 - k), ldbm,
                                             A(B->nt - 1 - k, B->nt - 1 - n), ldak,
                                             lalpha, B(m, B->nt - 1 - n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, B(m, B->nt - 1 - k));
                    }
                    for (n = k + 1; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, A(B->nt - 1 - k, B->nt - 1 - n));
                    }
                }
            }
                /*
                 *  HicmaRight / HicmaLower / Hicma[Conj]Trans
                 */
            else {
                for (k = 0; k < B->nt; k++) {
                    tempkn = k == B->nt - 1 ? B->n - k * B->nb : B->nb;
                    ldak = BLKLDD(A, k);
                    for (m = 0; m < B->mt; m++) {
                        tempmm = m == B->mt - 1 ? B->m - m * B->mb : B->mb;
                        ldbm = BLKLDD(B, m);
//                        HICMA_TASK_hcore_dtrsm(
//                                &options,
//                                side, uplo, trans, diag,
//                                tempmm, tempkn, A->mb,
//                                alpha, A(k, k), ldak,  /* lda * tempkn */
//                                B(m, k), ldbm); /* ldb * tempkn */

                        HICMA_TASK_dtrsm(&options, side, uplo, trans, diag,
                                               tempmm, tempkn, alpha, A(k, k), ldak,
                                               B(m, k), ldbm);

                        HICMA_RUNTIME_data_flush(sequence, A(k, k));

                        for (n = k + 1; n < B->nt; n++) {
                            tempnn = n == B->nt - 1 ? B->n - n * B->nb : B->nb;
                            ldan = BLKLDD(A, n);
//                            HICMA_TASK_dgemm(
//                                    &options,
//                                    HicmaNoTrans, trans,
//                                    tempmm, tempnn, B->mb, A->mb,
//                                    minvalpha, B(m, k), ldbm,  /* ldb  * tempkn */
//                                    A(n, k), ldan, /* ldan * tempkn */
//                                    done, B(m, n), ldbm); /* ldb  * tempnn */

                            HICMA_TASK_dgemm(&options, HicmaNoTrans, trans,
                                             tempmm, tempnn, B->mb, minvalpha,
                                             B(m, k), ldbm,
                                             A(n, k), ldan,
                                             done, B(m, n), ldbm);

                        }
                        HICMA_RUNTIME_data_flush(sequence, B(m, k));
                    }
                    for (n = k + 1; n < B->nt; n++) {
                        HICMA_RUNTIME_data_flush(sequence, A(n, k));
                    }

                }
            }
        }
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicma);
}

