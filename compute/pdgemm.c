/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 *
 * @file pdgemm.c
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Ali Charara
 * @author Kadir Akbudak
 * @date 2019-10-21
 *
 **/

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 **/

/**
 *
 * @file pdgemm.c
 *
 *  MORSE auxiliary routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 2.5.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2019-10-21
 *
 **/

#include "hicma_common.h"
#include <hicma.h>
#include <control/common.h>
#include "hicma_runtime_d.h"
#include "coreblas/hicma_lapacke.h"


#include "control/hicma_config.h"

// #define SYNCHRONOUS

#define A(m, n) AUV,  m,  n
#define B(m, n) BUV,  m,  n
#define C(m, n) CUV,  m,  n

#define AUV(m, n) AUV, Ark,  m,  n
#define BUV(m, n) BUV, Brk,  m,  n
#define CUV(m, n) CUV, Crk,  m,  n

#include "hicma.h"

int reorder_inner_products =0;
//extern idrank ***iporder; //inner product order
idrank ***iporder = NULL;
//int reorder_inner_products =0;
/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 **/
void hicma_pdgemm(HICMA_enum transA, HICMA_enum transB,
                  double alpha, HICMA_desc_t *AUV, HICMA_desc_t *Ark,
        // HICMA_Complex64_t alpha, HICMA_desc_t *AUV, HICMA_desc_t *Ark,
                  HICMA_desc_t *BUV, HICMA_desc_t *Brk,
                  double beta, HICMA_desc_t *CUV, HICMA_desc_t *Crk,
        // HICMA_Complex64_t beta,  HICMA_desc_t *CUV, HICMA_desc_t *Crk,
                  HICMA_sequence_t *sequence, HICMA_request_t *request,
                  int rk, int maxrk, double acc) {
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;
    size_t ws_host = 0;
    size_t ws_worker = 0;


    double dbeta;
    double done = (double) 1.0;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);

    ws_worker =  //FIXME tentative size. FInd exact size. I think syrk uses less memory
            //Ali says: this workspace need to be fixed, not all tasks below need it nor need that much
            2 * CUV->mb * 2 * maxrk //CUV clone
            + 2 * CUV->mb       // qrtauA qrtauB
            + maxrk * maxrk    // qrb_aubut  AcolBcolT
            + 2 * CUV->mb * maxrk // newU newV
            + (2 * maxrk) * (2 * maxrk)    // svd_rA     _rA
            //+ maxrk * maxrk    // svd_rB     _rB   I assume that use_trmm=1 so I commented out
            //+ maxrk * maxrk    // svd_T      _T    I assume that use_trmm=1 so I commented out
            + (2 * maxrk)          // sigma
#ifdef HCORE_GEMM_USE_ORGQR
        + CUV->mb * 2*maxrk // newUV gemms
#endif
            ;
    if (HICMA_get_use_fast_hcore_gemm() == 1) {
        double work_query;
        int lwork = -1;
        int info = LAPACKE_dgesvd_work(LAPACK_COL_MAJOR, 'A', 'A',
                                       2 * maxrk, 2 * maxrk,
                                       NULL, 2 * maxrk,
                                       NULL,
                                       NULL, 2 * maxrk,
                                       NULL, 2 * maxrk, &work_query, lwork);
        lwork = (int) work_query;
        ws_worker += lwork; // superb
    } else {
        ws_worker += (2 * maxrk); // superb
    }

    ws_worker *= sizeof(double);
    HICMA_RUNTIME_options_ws_alloc(&options, ws_worker, ws_host);

    for (m = 0; m < CUV->mt; m++) {
        tempmm = m == CUV->mt - 1 ? CUV->m - m * CUV->mb : CUV->mb;
        ldcm = BLKLDD(CUV, m);
        for (n = 0; n < CUV->nt; n++) {
            tempnn = n == CUV->nt - 1 ? CUV->n - n * CUV->nb : CUV->nb;
            /*
             *  A: HicmaNoTrans / B: HicmaNoTrans
             */
            if (transA == HicmaNoTrans) {
                ldam = BLKLDD(AUV, m);
                if (transB == HicmaNoTrans) {
                    for (int _k = 0; _k < AUV->nt; _k++) {
                        if (reorder_inner_products == 1) {
                            k = iporder[m][n][_k].id;
                        } else {
                            k = _k;
                        }
                        tempkn = k == AUV->nt - 1 ? AUV->n - k * AUV->nb : AUV->nb;
                        ldbk = BLKLDD(BUV, k);
                        dbeta = k == 0 ? beta : done;
                        HICMA_TASK_hcore_dgemm(
                                &options,
                                transA, transB,
                                tempmm, tempnn,
                                alpha, AUV(m, k), ldam,  /* lda * Z */
                                BUV(k, n), ldbk,  /* ldb * Y */
                                dbeta, CUV(m, n), ldcm, /* ldc * Y */
                                rk, maxrk, acc);
#ifdef SYNCHRONOUS
                        HICMA_RUNTIME_barrier(hicma);
#endif
                    }
                    //HICMA_RUNTIME_barrier(hicma); //for forcing the execution order to be same as submission order
                }
                    /*
                     *  A: HicmaNoTrans / B: Hicma[Conj]Trans
                     */
                else {
                    ldbn = BLKLDD(BUV, n);
                    for (k = 0; k < AUV->nt; k++) {
                        tempkn = k == AUV->nt - 1 ? AUV->n - k * AUV->nb : AUV->nb;
                        dbeta = k == 0 ? beta : done;
                        HICMA_TASK_hcore_dgemm(
                                &options,
                                transA, transB,
                                tempmm, tempnn,
                                alpha, AUV(m, k), ldam,  /* lda * Z */
                                BUV(n, k), ldbn,  /* ldb * Z */
                                dbeta, CUV(m, n), ldcm, /* ldc * Y */
                                rk, maxrk, acc);
#ifdef SYNCHRONOUS
                        HICMA_RUNTIME_barrier(hicma);
#endif
                    }
                }
            }
                /*
                 *  A: Hicma[Conj]Trans / B: HicmaNoTrans
                 */
            else {
                if (transB == HicmaNoTrans) {
                    for (k = 0; k < AUV->mt; k++) {
                        tempkm = k == AUV->mt - 1 ? AUV->m - k * AUV->mb : AUV->mb;
                        ldak = BLKLDD(AUV, k);
                        ldbk = BLKLDD(BUV, k);
                        dbeta = k == 0 ? beta : done;
                        HICMA_TASK_hcore_dgemm(
                                &options,
                                transA, transB,
                                tempmm, tempnn,
                                alpha, AUV(k, m), ldak,  /* lda * X */
                                BUV(k, n), ldbk,  /* ldb * Y */
                                dbeta, CUV(m, n), ldcm, /* ldc * Y */
                                rk, maxrk, acc);
#ifdef SYNCHRONOUS
                        HICMA_RUNTIME_barrier(hicma);
#endif
                    }
                }
                    /*
                     *  A: Hicma[Conj]Trans / B: Hicma[Conj]Trans
                     */
                else {
                    ldbn = BLKLDD(BUV, n);
                    for (k = 0; k < AUV->mt; k++) {
                        tempkm = k == AUV->mt - 1 ? AUV->m - k * AUV->mb : AUV->mb;
                        ldak = BLKLDD(AUV, k);
                        dbeta = k == 0 ? beta : done;
                        HICMA_TASK_hcore_dgemm(
                                &options,
                                transA, transB,
                                tempmm, tempnn,
                                alpha, AUV(k, m), ldak,  /* lda * X */
                                BUV(n, k), ldbn,  /* ldb * Z */
                                dbeta, CUV(m, n), ldcm, /* ldc * Y */
                                rk, maxrk, acc);
#ifdef SYNCHRONOUS
                        HICMA_RUNTIME_barrier(hicma);
#endif
                    }
                }
            }
            HICMA_RUNTIME_data_flush(sequence, C(m, n));
        }
        if (transA == HicmaNoTrans) {
            for (k = 0; k < AUV->nt; k++) {
                //HICMA_TASK_dataflush( &options, A(m, k) );
                HICMA_RUNTIME_data_flush(sequence, A(m, k));
            }
        } else {
            for (k = 0; k < AUV->mt; k++) {
                /*HICMA_TASK_dataflush( &options, A(k, m) );*/
                HICMA_RUNTIME_data_flush(sequence, A(k, m));
            }
        }
        /*for (n = 0; n < CUV->nt; n++) {*/
        /*HICMA_TASK_dataflush( &options, C(m, n) );*/
        /*}*/
    }
    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicma);
    //HICMA_TASK_dataflush_all(); removed in newer chameleon
}
