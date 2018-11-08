/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 *
 * @file pzgemm.c
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Ali Charara
 * @author Kadir Akbudak
 * @date 2018-11-08
 *
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 */

/*
 *
 * @file pzgemm.c
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
 * @date 2018-11-08
 * @precisions normal z -> s d c
 *
 **/
#include "hicma_common.h"
#include "morse.h"
#include "control/common.h"
#include "hicma_runtime_z.h"
#include "coreblas/lapacke.h"

#include "control/hicma_config.h"

// #define SYNCHRONOUS

#define A(m, n) AUV,  m,  n
#define B(m, n) BUV,  m,  n
#define C(m, n) CUV,  m,  n

#define AUV(m, n) AUV, Ark,  m,  n
#define BUV(m, n) BUV, Brk,  m,  n
#define CUV(m, n) CUV, Crk,  m,  n

#include "hicma.h"

/**
 *  Parallel tile matrix-matrix multiplication - dynamic scheduling
 **/
void hicma_pzgemm(MORSE_enum transA, MORSE_enum transB,
                  double alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
                  // MORSE_Complex64_t alpha, MORSE_desc_t *AUV, MORSE_desc_t *Ark,
                                           MORSE_desc_t *BUV, MORSE_desc_t *Brk,
                  double beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
                  // MORSE_Complex64_t beta,  MORSE_desc_t *CUV, MORSE_desc_t *Crk,
                  MORSE_sequence_t *sequence, MORSE_request_t *request,
                  int rk, int maxrk, double acc)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int m, n, k;
    int ldam, ldak, ldbn, ldbk, ldcm;
    int tempmm, tempnn, tempkn, tempkm;
    size_t ws_host   = 0;
    size_t ws_worker   = 0;


    double zbeta;
    double zone = (double)1.0;
    // MORSE_Complex64_t zbeta;
    // MORSE_Complex64_t zone = (MORSE_Complex64_t)1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);

    ws_worker =  //FIXME tentative size. FInd exact size. I think syrk uses less memory
    //Ali says: this workspace need to be fixed, not all tasks below need it nor need that much
        2 * CUV->mb * 2 * maxrk //CUV clone
        + 2 * CUV->mb       // qrtauA qrtauB
        + maxrk * maxrk    // qrb_aubut  AcolBcolT
        + 2 * CUV->mb * maxrk // newU newV
        + (2*maxrk) * (2*maxrk)    // svd_rA     _rA
        //+ maxrk * maxrk    // svd_rB     _rB   I assume that use_trmm=1 so I commented out
        //+ maxrk * maxrk    // svd_T      _T    I assume that use_trmm=1 so I commented out
        + (2*maxrk)          // sigma
        #ifdef HCORE_GEMM_USE_ORGQR
        + CUV->mb * 2*maxrk // newUV gemms
        #endif
        ;
    if(HICMA_get_use_fast_hcore_zgemm() == 1){
        double work_query;
        int lwork = -1;
        int info = LAPACKE_dgesvd_work( LAPACK_COL_MAJOR, 'A', 'A',
            2*maxrk, 2*maxrk,
            NULL, 2*maxrk,
            NULL,
            NULL, 2*maxrk,
            NULL, 2*maxrk, &work_query, lwork );
        lwork = (int)work_query;
        ws_worker += lwork; // superb
    }else{
        ws_worker += (2*maxrk); // superb
    }

    ws_worker *= sizeof(double); //FIXME use MORSE_Complex64_t
    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );

    for (m = 0; m < CUV->mt; m++) {
        tempmm = m == CUV->mt-1 ? CUV->m-m*CUV->mb : CUV->mb;
        ldcm = BLKLDD(CUV, m);
        for (n = 0; n < CUV->nt; n++) {
            tempnn = n == CUV->nt-1 ? CUV->n-n*CUV->nb : CUV->nb;
            /*
             *  A: MorseNoTrans / B: MorseNoTrans
             */
            if (transA == MorseNoTrans) {
                ldam = BLKLDD(AUV, m);
                if (transB == MorseNoTrans) {
                    for (k = 0; k < AUV->nt; k++) {
                        tempkn = k == AUV->nt-1 ? AUV->n-k*AUV->nb : AUV->nb;
                        ldbk = BLKLDD(BUV, k);
                        zbeta = k == 0 ? beta : zone;
                        HICMA_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, 
                            alpha, AUV(m, k), ldam,  /* lda * Z */
                                   BUV(k, n), ldbk,  /* ldb * Y */
                            zbeta, CUV(m, n), ldcm, /* ldc * Y */
                            rk, maxrk, acc);
				                #ifdef SYNCHRONOUS
				                RUNTIME_barrier(morse);
				                #endif
                    }
                }
                /*
                 *  A: MorseNoTrans / B: Morse[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(BUV, n);
                    for (k = 0; k < AUV->nt; k++) {
                        tempkn = k == AUV->nt-1 ? AUV->n-k*AUV->nb : AUV->nb;
                        zbeta = k == 0 ? beta : zone;
                        HICMA_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, 
                            alpha, AUV(m, k), ldam,  /* lda * Z */
                                   BUV(n, k), ldbn,  /* ldb * Z */
                            zbeta, CUV(m, n), ldcm, /* ldc * Y */
                            rk, maxrk, acc);
		                #ifdef SYNCHRONOUS
		                RUNTIME_barrier(morse);
		                #endif
                    }
                }
            }
            /*
             *  A: Morse[Conj]Trans / B: MorseNoTrans
             */
            else {
                if (transB == MorseNoTrans) {
                    for (k = 0; k < AUV->mt; k++) {
                        tempkm = k == AUV->mt-1 ? AUV->m-k*AUV->mb : AUV->mb;
                        ldak = BLKLDD(AUV, k);
                        ldbk = BLKLDD(BUV, k);
                        zbeta = k == 0 ? beta : zone;
                        HICMA_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, 
                            alpha, AUV(k, m), ldak,  /* lda * X */
                                   BUV(k, n), ldbk,  /* ldb * Y */
                            zbeta, CUV(m, n), ldcm, /* ldc * Y */
                            rk, maxrk, acc);
		                #ifdef SYNCHRONOUS
		                RUNTIME_barrier(morse);
		                #endif
                    }
                }
                /*
                 *  A: Morse[Conj]Trans / B: Morse[Conj]Trans
                 */
                else {
                    ldbn = BLKLDD(BUV, n);
                    for (k = 0; k < AUV->mt; k++) {
                        tempkm = k == AUV->mt-1 ? AUV->m-k*AUV->mb : AUV->mb;
                        ldak = BLKLDD(AUV, k);
                        zbeta = k == 0 ? beta : zone;
                        HICMA_TASK_zgemm(
                            &options,
                            transA, transB,
                            tempmm, tempnn, 
                            alpha, AUV(k, m), ldak,  /* lda * X */
                                   BUV(n, k), ldbn,  /* ldb * Z */
                            zbeta, CUV(m, n), ldcm, /* ldc * Y */
                            rk, maxrk, acc);
		                #ifdef SYNCHRONOUS
		                RUNTIME_barrier(morse);
		                #endif
                    }
                }
            }
            RUNTIME_data_flush( sequence, C(m, n) );
        }
        if (transA == MorseNoTrans) {
            for (k = 0; k < AUV->nt; k++) {
                //MORSE_TASK_dataflush( &options, A(m, k) );
                RUNTIME_data_flush( sequence, A(m, k) );
            }
        } else {
            for (k = 0; k < AUV->mt; k++) {
                /*MORSE_TASK_dataflush( &options, A(k, m) );*/
                RUNTIME_data_flush( sequence, A(k, m) );
            }
        }
        /*for (n = 0; n < CUV->nt; n++) {*/
            /*MORSE_TASK_dataflush( &options, C(m, n) );*/
        /*}*/
    }
    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
    //MORSE_TASK_dataflush_all(); removed in newer chameleon 
}
