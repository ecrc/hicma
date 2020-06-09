/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file pzpotrf.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */
/*
 *
 * file pzpotrf.c
 *
 * MORSE auxiliary routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * author Jakub Kurzak
 * author Hatem Ltaief
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * author Florent Pruvost
 * date 2010-11-15
 *
 **/
#include "morse.h"
#include "hicma.h"
#include "hicma_common.h"
#include "control/common.h"
#include "hicma_runtime_z.h"
#include "coreblas/lapacke.h"

#include "control/hicma_config.h"
#include <stdio.h>

extern int store_only_diagonal_tiles;
extern int print_index;
int pzgetrf_print_index = 0;
extern int print_mat;
extern int run_org;
int extra_barrier = 0;
/***************************************************************************//**
 *  Parallel tile Cholesky factorization - dynamic scheduling
 **/
void hicma_pzgetrf(MORSE_enum uplo,
        MORSE_desc_t *AUV,
        MORSE_desc_t *AD,
        MORSE_desc_t *Ark,
        MORSE_sequence_t *sequence, MORSE_request_t *request,
        int rk, int maxrk, double acc)
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int k, m, n, i, j;
    size_t ws_host   = 0;
    size_t ws_worker   = 0;

    MORSE_Complex64_t zone  = (MORSE_Complex64_t)1.0;
    MORSE_Complex64_t mzone = (MORSE_Complex64_t)-1.0;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;
    RUNTIME_options_init(&options, morse, sequence, request);


/*#ifdef CHAMELEON_USE_MAGMA*/
    /*if (0) [> Disable the workspace as long as it is is not used (See StarPU codelet) <]*/
    /*{*/
        /*int nb = MORSE_IB; [> Approximate nb for simulation <]*/
/*#if !defined(CHAMELEON_SIMULATION)*/
        /*nb = magma_get_zpotrf_nb(AD->nb);*/
/*#endif*/
        /*ws_host = sizeof(double)*nb*nb;*/
    /*}*/
/*#endif*/
    //printf("%s %s %d maxrank=%d\n", __FILE__, __func__, __LINE__, maxrk);
    ws_worker =  //FIXME tentative size. FInd exact size. I think syrk uses less memory
        //Ali says: this workspace need to be fixed, not all tasks below need it nor need that much
        2 * AD->mb * 2 * maxrk   // for copying CU and CV into temporary buffer instead of using CUV itself. There is 2*maxrk because these buffers will be used to put two U's side by side
        + 2 * AD->mb       // qrtauA qrtauB
        + maxrk * maxrk    // qrb_aubut  AcolBcolT
        + 2 * chameleon_max(AD->mb, maxrk) * maxrk // newU newV
        + (2*maxrk) * (2*maxrk)    // svd_rA     _rA
        + maxrk * maxrk    // svd_rB     _rB   TRMM will NOT be used so this sentence is invalid:I assume that use_trmm=1 so I commented out
        + maxrk * maxrk    // svd_T      _T    TRMM will NOT be used so this sentence is invalid:I assume that use_trmm=1 so I commented out
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

    ws_worker *= sizeof(MORSE_Complex64_t); 
    //FIXME add ws_worker and ws_host calculation from compute/pzgeqrf.c when GPU/MAGMA is supported
    RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );


    /*
     *  MorseLower
     */
    if (uplo == MorseLower) {
        for (k = 0; k < AD->mt; k++) {
            RUNTIME_iteration_push(morse, k);

            int tempkmd = k == AD->mt-1 ? AD->m-k*AD->mb : AD->mb;
            int ldakd = BLKLDD(AD, k);

            //options.priority = 2*AD->mt - 2*k;
            options.priority = 5;
            if(pzgetrf_print_index){
                printf("GETRF\t|panel/k:%d tempkmd:%d ldakd:%d\n", k, tempkmd, ldakd);
            }
            int ADicol;
            if(store_only_diagonal_tiles == 1){
                ADicol = 0;
            } else {
                ADicol = k;
            }
           if(1)HICMA_TASK_zgetrf(
                   &options,
                   tempkmd, AD->mb,
                   AD, k, ADicol, ldakd, 0);
                  
/*            HICMA_TASK_zgetrf(
                    &options,
                    tempkmd, AD->mb,
                    AD, k, ADicol, ldakd, 0);
*/
            for (i = k+1; i < AD->mt; i++) {
                int ldamuv = BLKLDD(AUV, i);

                //options.priority = 2*AD->mt - 2*k - m;
                options.priority = 4;
                if(pzgetrf_print_index){
                    printf("TRSM\t|m:%d k:%d ldakd:%d ldamuv:%d\n", i, k, ldakd, ldamuv);
                }
                /*
                 * X D^t = U V^t
                 u X = U V^t * inv(D^t)
                 * X = U (inv(D) * V )^t
                 * X = U * (trsm (lower, left, notranspose, D, V) )^t
                 * X = U * newV^t
                 */
                if(1)HICMA_TASK_ztrsm(
                        &options,
                        MorseLeft, MorseUpper, MorseTrans, MorseNonUnit,
                        tempkmd, // number of rows of the diagonal block  
                        zone, AD, k, ADicol, 
                        ldakd,
                        AUV, i, k, 
                        ldamuv, 
                        Ark);
            }

              /* HICMA_TASK_ztrsmu(
                        &options,
                        MorseLeft, MorseLower, MorseNoTrans, MorseUnit,
                        tempkmd, // number of rows of the diagonal block  
                        zone, AD, k, ADicol, 
                        ldakd,
                        AUV, m, k, 
                        ldamuv, 
                        Ark); */

            //MORSE_TASK_dataflush( &options, AV, k, k );
            //MORSE_TASK_dataflush( &options, AUV, k, k );
            //RUNTIME_data_flush( sequence, AUV, k, k); 

            for (j = k+1; j < AD->mt; j++) {
                int ldamuv = BLKLDD(AUV, k);

                //options.priority = 2*AD->mt - 2*k - m;
                options.priority = 3;
                if(pzgetrf_print_index){
                    printf("TRSMU\t|k:%d j:%d ldakd:%d ldamuv:%d\n", k, j, ldakd, ldamuv);
                }
                /*
                 * X D^t = U V^t
                 * X = U V^t * inv(D^t)
                 * X = U (inv(D) * V )^t
                 * X = U * (trsm (lower, left, notranspose, D, V) )^t
                 * X = U * newV^t
                 */
                if(1)HICMA_TASK_ztrsmu(
                        &options,
                        MorseLeft, MorseLower, MorseNoTrans, MorseUnit,
                        tempkmd, // number of rows of the diagonal block
                        zone, AD, k, ADicol,
                        ldakd,
                        AUV, k,j,
                        ldamuv,
                        Ark);
            }
            //RUNTIME_data_flush( sequence, AUV, k, k);
           RUNTIME_data_flush( sequence, AD, k, ADicol);
            for (i = k+1; i < AD->mt; i++) {
                int tempnnd = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
                int ldand = BLKLDD(AD, i);
                int ldanuv = BLKLDD(AUV, i);
                //int ldanuv = BLKLDD(BUV, i;
                int ADicol;
                if(store_only_diagonal_tiles == 1){
                    ADicol = 0;
                } else {
                    ADicol = i;
                }

                //options.priority = 2*AD->mt - 2*k - n;
                options.priority = 2;
                if(1)HICMA_TASK_zgemm_cd(
                        &options,
                        tempnnd,   0,
                        -1.0,
                        AUV, ldanuv,
                        Ark,
                        i, k,
                        AUV, ldanuv,
                        Ark,
                        k,i,
                        1.0,
                        AD, ldand,
                        i, ADicol
                        );

                for (j =i+1; j < AD->mt; j++) {
                    int tempmmuv = j == AUV->mt-1 ? AUV->m - j*AUV->mb : AUV->mb;
                    int ldamuv = BLKLDD(AUV, j);

                    //options.priority = 2*AD->mt - 2*k - n - m;
                    options.priority =1; // 2*AD->mt - 2*i - j;
                    if(pzgetrf_print_index ){
                        printf("GEMM\t|A(%d,%d)=A(%d,%d)-A(%d,%d)*A(%d,%d) ldamuv:%d tempmmuv:%d\n", 
                              //m,n,m,n,m,k,n,k,
                              j,i,j,i,j,k,i,k,
                              ldamuv, tempmmuv);
                        printf("GEMM-1\t|C(%d,%d)=C(%d,%d)-A(%d,%d)*B(%d,%d)\n",i,j,i,j,i,k,k,j);
                        printf("GEMM-2\t|C(%d,%d)=C(%d,%d)-A(%d,%d)*B(%d,%d)\n",j,i,j,i,j,k,k,i);
                    }
                    if(1) HICMA_TASK_zgemm(
                            &options,
                            MorseNoTrans, MorseNoTrans,
                            tempmmuv,
                            tempmmuv,
                            mzone, AUV, Ark, i, k, ldamuv,
                            AUV, Ark, k, j, ldamuv,
                            zone, AUV, Ark, i, j, ldamuv,
                            rk, maxrk, acc);

                    if(1) HICMA_TASK_zgemm(
                            &options,
                            MorseNoTrans, MorseNoTrans,
                            tempmmuv,
                            tempmmuv,
                            mzone, AUV, Ark, j, k, ldamuv,
                            AUV, Ark, k, i, ldamuv,
                            zone, AUV, Ark, j, i, ldamuv,
                            rk, maxrk, acc);
                }               
                //MORSE_TASK_dataflush( &options, AUV, n, k );
                //RUNTIME_data_flush( sequence, AUV, i, k);
                RUNTIME_data_flush( sequence, AUV, i, k);
                RUNTIME_data_flush( sequence, AUV, k, i);
                RUNTIME_data_flush( sequence, Ark, i, k);
                RUNTIME_data_flush( sequence, Ark, k, i);
            }
            RUNTIME_iteration_pop(morse);

            if(extra_barrier){
//                RUNTIME_barrier(morse);
            }
        }
    }
    /*
     *  MorseUpper
     */


    RUNTIME_options_ws_free(&options);
    RUNTIME_options_finalize(&options, morse);
}

