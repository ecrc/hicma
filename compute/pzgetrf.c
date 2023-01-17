/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file pzgetrf.c
 *
 *  HiCMA auxiliary routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
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
 * file pzgetrf.c
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

#include <stdio.h>

#include <hicma.h>
#include <control/common.h>
#include <include/hicma_runtime_z.h>

#include <coreblas/hicma_lapacke.h>

extern int store_only_diagonal_tiles;
extern int print_index;
int pzgetrf_print_index = 0;
extern int print_mat;
extern int run_org;
/***************************************************************************//**
 *  Parallel tile LU factorization - dynamic scheduling
 **/
void hicma_pzgetrf(HICMA_enum uplo,
        HICMA_desc_t *AUV,
        HICMA_desc_t *AD,
        HICMA_desc_t *Ark,
        HICMA_sequence_t *sequence, HICMA_request_t *request,
        int rk, int maxrk, double acc)
{
    HICMA_context_t *hicma;
    HICMA_option_t options;

    int k, m, n, i, j;
    size_t ws_host   = 0;
    size_t ws_worker   = 0;

    HICMA_Complex64_t zone  = (HICMA_Complex64_t)1.0;
    HICMA_Complex64_t mzone = (HICMA_Complex64_t)-1.0;

    hicma = hicma_context_self();
    if (sequence->status != HICMA_SUCCESS)
        return;
    HICMA_RUNTIME_options_init(&options, hicma, sequence, request);


/*#ifdef CHAMELEON_USE_MAGMA*/
    /*if (0) [> Disable the workspace as long as it is is not used (See StarPU codelet) <]*/
    /*{*/
        /*int nb = HICMA_IB; [> Approximate nb for simulation <]*/
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
        + 2 * hicma_max(AD->mb, maxrk) * maxrk // newU newV
        + (2*maxrk) * (2*maxrk)    // svd_rA     _rA
        + maxrk * maxrk    // svd_rB     _rB   TRMM will NOT be used so this sentence is invalid:I assume that use_trmm=1 so I commented out
        + maxrk * maxrk    // svd_T      _T    TRMM will NOT be used so this sentence is invalid:I assume that use_trmm=1 so I commented out
        + (2*maxrk)          // sigma
        #ifdef HCORE_GEMM_USE_ORGQR
        + CUV->mb * 2*maxrk // newUV gemms
        #endif
        ;
    if(HICMA_get_use_fast_hcore_gemm() == 1){
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

    ws_worker *= sizeof(HICMA_Complex64_t);
    //FIXME add ws_worker and ws_host calculation from compute/pzgeqrf.c when GPU/MAGMA is supported
    HICMA_RUNTIME_options_ws_alloc( &options, ws_worker, ws_host );


    for (k = 0; k < AD->mt; k++) {
        HICMA_RUNTIME_iteration_push(hicma, k);

        int tempkmd = k == AD->mt-1 ? AD->m-k*AD->mb : AD->mb;
        int ldakd = BLKLDD(AD, k);

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
        HICMA_TASK_zgetrf(
                &options,
                tempkmd, AD->mb,
                AD, k, ADicol, ldakd, 0);
        
        for (i = k+1; i < AD->mt; i++) {
            int ldamuv = BLKLDD(AUV, i);

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
            HICMA_TASK_hcore_ztrsm(
                    &options,
                    HicmaLeft, HicmaUpper, HicmaTrans, HicmaNonUnit,
                    tempkmd, // number of rows of the diagonal block  
                    zone, AD, k, ADicol, 
                    ldakd,
                    AUV, i, k, 
                    ldamuv, 
                    Ark);
        }

        for (j = k+1; j < AD->mt; j++) {
            int ldamuv = BLKLDD(AUV, k);

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
            HICMA_TASK_hcore_ztrsmu(
                    &options,
                    HicmaLeft, HicmaLower, HicmaNoTrans, HicmaUnit,
                    tempkmd, // number of rows of the diagonal block
                    zone, AD, k, ADicol,
                    ldakd,
                    AUV, k,j,
                    ldamuv,
                    Ark);
        }

        HICMA_RUNTIME_data_flush( sequence, AD, k, ADicol);
        for (i = k+1; i < AD->mt; i++) {
            int tempnnd = i == AD->mt-1 ? AD->m-i*AD->mb : AD->mb;
            int ldand = BLKLDD(AD, i);
            int ldanuv = BLKLDD(AUV, i);
            int ADicol;
            if(store_only_diagonal_tiles == 1){
                ADicol = 0;
            } else {
                ADicol = i;
            }

            options.priority = 2;
            HICMA_TASK_zgemm_cd(
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

                options.priority =1; // 2*AD->mt - 2*i - j;
                if(pzgetrf_print_index ){
                    printf("GEMM\t|A(%d,%d)=A(%d,%d)-A(%d,%d)*A(%d,%d) ldamuv:%d tempmmuv:%d\n", 
                            //m,n,m,n,m,k,n,k,
                            j,i,j,i,j,k,i,k,
                            ldamuv, tempmmuv);
                    printf("GEMM-1\t|C(%d,%d)=C(%d,%d)-A(%d,%d)*B(%d,%d)\n",i,j,i,j,i,k,k,j);
                    printf("GEMM-2\t|C(%d,%d)=C(%d,%d)-A(%d,%d)*B(%d,%d)\n",j,i,j,i,j,k,k,i);
                }
                HICMA_TASK_hcore_zgemm(
                        &options,
                        HicmaNoTrans, HicmaNoTrans,
                        tempmmuv,
                        tempmmuv,
                        mzone, AUV, Ark, i, k, ldamuv,
                        AUV, Ark, k, j, ldamuv,
                        zone, AUV, Ark, i, j, ldamuv,
                        rk, maxrk, acc);

                HICMA_TASK_hcore_zgemm(
                        &options,
                        HicmaNoTrans, HicmaNoTrans,
                        tempmmuv,
                        tempmmuv,
                        mzone, AUV, Ark, j, k, ldamuv,
                        AUV, Ark, k, i, ldamuv,
                        zone, AUV, Ark, j, i, ldamuv,
                        rk, maxrk, acc);
            }               
            HICMA_RUNTIME_data_flush( sequence, AUV, i, k);
            HICMA_RUNTIME_data_flush( sequence, AUV, k, i);
            HICMA_RUNTIME_data_flush( sequence, Ark, i, k);
            HICMA_RUNTIME_data_flush( sequence, Ark, k, i);
        }
        HICMA_RUNTIME_iteration_pop(hicma);
    }

    HICMA_RUNTIME_options_ws_free(&options);
    HICMA_RUNTIME_options_finalize(&options, hicma);
}

