/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zgemm.c
 *
 *  HiCMA HCORE kernels
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/
#include "coreblas/include/coreblas.h"
#include "coreblas/lapacke.h"
#include <assert.h>
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif

#include "control/hicma_config.h"

//FIXME PREVIOUS DECLARION OF CBLAS_SADDR ~/hicma-dev/chameleon/build/include/chameleon/coreblas/include/coreblas.h
#undef CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

// #define DBG_MSG

#ifdef DBG_MSG
#define ECHO_I(_val) printf("%s(%d) ", #_val, (_val));
#define ECHO_f(_val) printf("%s(%e) ", #_val, (_val));
#define ECHO_LN printf("\n");
#else
#define ECHO_I(_val)
#define ECHO_f(_val)
#define ECHO_LN
#endif

#ifdef HCORE_GEMM_USE_KBLAS_ACA
extern int kblas_ACAf( int m, int n,
          double* A, int lda,
          double* U, int ldu,
          double* V, int ldv,
          double* S,
          double maxacc, int maxrk,
          double* acc, int* rk);
#endif

extern int use_trmm;
extern int use_scratch;
extern int gemm_print_index;
extern int gemm_print_mat;
extern void hc_printmat(double * A, int m, int n, int ld);
/***************************************************************************//*
 *
 * @ingroup CORE_double
 *
 **/

void HCORE_zgemm_fast(MORSE_enum transA, int transB,
        int M, int N,
        double alpha,
        double *AU,
        double *AV,
        double *Ark,
        int LDA,
        double *BU,
        double *BV,
        double *Brk,
        int LDB,
        double beta,
        double *CU,
        double *CV,
        double *Crk,
        int LDC,
        int rk,
        int maxrk,
        double acc,
        double* d_work
        )
{do{
    // printf("%d GEMM %p\n", MORSE_My_Mpi_Rank(), d_work);
    assert(use_trmm == 1);
    assert(use_scratch == 1);

    int ws_needed = 0;

    // cudaStat = cudaStreamSynchronize( cuda_stream );
    // assert(cudaSuccess == cudaStat);

    int _Ark = (int)(Ark[0]); ECHO_I(_Ark);
    int _Brk = (int)(Brk[0]); ECHO_I(_Brk);
    int _Crk = (int)(Crk[0]); ECHO_I(_Crk);
    int old_Crk = _Crk;
    // if(gemm_print_index) printf("%d, _Ark %d, _Brk %d, _Crk %d\n", __LINE__, _Ark, _Brk, _Crk);

    int _M = M; int _N = N;  ECHO_I(_M);
    double* _CU = CU; int ld_CU = LDC; ECHO_I(ld_CU);
    double* _CV = CV; int ld_CV = LDC; ECHO_I(ld_CV);
    double* _AU = AU; int ld_AU = LDA; ECHO_I(ld_AU);
    double* _AV = AV; int ld_AV = LDA; ECHO_I(ld_AV);
    double* _BU = BU; int ld_BU = LDB; ECHO_I(ld_BU);
    double* _BV = BV; int ld_BV = LDB; ECHO_I(ld_BV);
    int rank = rk; ECHO_I(rank); ECHO_I(maxrk)

    char chall = 'A';

    int use_CUV_clone = 0;
    double* _CU_save = _CU;
    double* _CV_save = _CV;
    int ld_CU_save = ld_CU;
    int ld_CV_save = ld_CV;

    int CUV_ncols = _Crk + _Ark; ECHO_I(CUV_ncols);

    if((CUV_ncols > maxrk)){
        double* CUclone = NULL;
        int ld_CUclone = _M;
        double* CVclone = NULL;
        int ld_CVclone = _M;
        size_t CUclone_nelm =  _M * 2 * maxrk;
        size_t CVclone_nelm =  _M * 2 * maxrk;

        use_CUV_clone = 1;
        CUclone = d_work;
        d_work += CUclone_nelm;
        ws_needed += CUclone_nelm;
        CVclone = d_work;
        d_work += CVclone_nelm;
        ws_needed += CVclone_nelm;
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CU, &ld_CU,
                CUclone, &ld_CUclone);
        LAPACK_dlacpy(&chall,
                &_M, &_Crk,
                _CV, &ld_CV,
                CVclone, &ld_CVclone);
        _CU = CUclone;
        _CV = CVclone;
        ld_CU = ld_CUclone;
        ld_CV = ld_CVclone;
    }
    int incOne = 1;
    double d_one = 1.0;
    double d_zero = 0.0;

    // TODO remove the abundant assumptions on matrices sizes and leading dimensions

    //=======================================================================================================
    // QR A

        //concat CU+AU
        int nelm_AU = _M * _Ark;  ECHO_I(nelm_AU); ECHO_I(_Crk*ld_CU);
        cblas_dcopy(nelm_AU, _AU, incOne,  &_CU[_Crk*ld_CU], incOne);

        if(alpha != d_one){
            cblas_dscal(nelm_AU, CBLAS_SADDR(alpha), &_CU[_Crk*ld_CU], incOne);
        }
        if(beta != d_one){
            ECHO_I(_M * _Crk);
            cblas_dscal(_M * _Crk, CBLAS_SADDR(beta), _CU, incOne);
        }

        double *qrtauA = d_work;
        size_t qrtauA_nelm = _M;
        d_work += qrtauA_nelm;
        ws_needed += qrtauA_nelm;
        assert(qrtauA != NULL);

        int info = LAPACKE_dgeqrf( LAPACK_COL_MAJOR, _M, CUV_ncols, _CU, ld_CU, qrtauA);

    //=======================================================================================================
    // QR B
        double* qrb_avtbv = d_work;
        size_t qrb_avtbv_nelm = maxrk * maxrk; ECHO_I(qrb_avtbv_nelm);
        d_work += qrb_avtbv_nelm;
        ws_needed += qrb_avtbv_nelm;

        //P = AV^T * BV
        cblas_dgemm(CblasColMajor,
                    CblasTrans, CblasNoTrans,
                    _Ark, _Brk, _M,
                    CBLAS_SADDR(d_one), _AV, ld_AV,
                                        _BV, ld_BV,
                    CBLAS_SADDR(d_zero), qrb_avtbv, maxrk);

        //G = P * BU^T <=> G^T = BU * P^T
        //CV = CV | G^T
        cblas_dgemm(CblasColMajor,
                    CblasNoTrans, CblasTrans,
                    _M, _Ark, _Brk,
                    CBLAS_SADDR(d_one), _BU, ld_BU,
                                        qrb_avtbv, maxrk,
                    CBLAS_SADDR(d_zero), &_CV[_Crk*ld_CV], ld_CV);

        // double* qrtauB = qrtauA + qrtauA_nelm;
        double* qrtauB = d_work;
        size_t  qrtauB_nelm = _M;
        d_work += qrtauB_nelm;
        ws_needed += qrtauB_nelm;
        assert(qrtauB != NULL);

        info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, _M, CUV_ncols, _CV, ld_CV, qrtauB);

    //=======================================================================================================
    //SVD
        double* rA = d_work;
        size_t rA_nelm = CUV_ncols * CUV_ncols; ECHO_I(rA_nelm);
        int ld_rA  = CUV_ncols;
        d_work += rA_nelm;
        ws_needed += rA_nelm;

        double* rB = _CV;
        int ld_rB  = ld_CV;
        // size_t rB_nelm = CV_ncols * CV_ncols;

        //copy rA from CU
        char chlow = 'L';
        LAPACK_dlaset(&chlow, &CUV_ncols, &CUV_ncols, &d_zero, &d_zero, rA, &ld_rA);
        char chup = 'U';
        LAPACK_dlacpy(&chup, &CUV_ncols, &CUV_ncols,
                        _CU, &ld_CU,
                        rA, &ld_rA);

        // rA = rA * rB^T
        cblas_dtrmm(CblasColMajor, CblasRight, CblasUpper, CblasTrans, CblasNonUnit,
                    CUV_ncols, CUV_ncols,
                    d_one, rB, ld_rB,
                           rA, ld_rA);

        int finalrank = -1, size_sigma = CUV_ncols; ECHO_I(size_sigma)
        double relacc = (acc);

        double* _T = rA;
        int ld_T = ld_rA;
        double* TU = d_work;
        #ifdef HCORE_GEMM_USE_ORGQR
            size_t TU_nelm = CUV_ncols * CUV_ncols; ECHO_I(TU_nelm)
            int ld_TU = CUV_ncols;
        #else
            size_t TU_nelm = _M * CUV_ncols; ECHO_I(TU_nelm)
            int ld_TU = _M;
        #endif
        d_work += TU_nelm;
        ws_needed += TU_nelm;

        double* TV = d_work;
        #ifdef HCORE_GEMM_USE_ORGQR
            size_t TV_nelm = CUV_ncols * CUV_ncols;
            int ld_TV = CUV_ncols;
        #else
            size_t TV_nelm = _M * CUV_ncols;
            #ifdef HCORE_GEMM_USE_KBLAS_ACA
                int ld_TV = _M;
            #else
                int ld_TV = CUV_ncols;
            #endif
        #endif
        ECHO_I(TV_nelm)
        d_work += TV_nelm;
        ws_needed += TV_nelm;

        double *d_sigma  = d_work;
        size_t d_sigma_nelm  = CUV_ncols; ECHO_I(d_sigma_nelm)
        d_work += d_sigma_nelm;
        ws_needed += d_sigma_nelm;

        #if defined HCORE_GEMM_USE_KBLAS_ACA

            double finalacc;
            kblas_ACAf( CUV_ncols, CUV_ncols,
                        _T, ld_T,
                        TU, ld_TU,
                        TV, ld_TV,
                        d_sigma,
                        relacc, rank,
                        &finalacc, &finalrank);
            ECHO_I(finalrank)
        #else
            double* svdsuperb = d_work;
            double work_query;
            int lwork = -1;
            info = LAPACKE_dgesvd_work( LAPACK_COL_MAJOR, 'A', 'A',
                                        CUV_ncols, CUV_ncols,
                                        NULL, CUV_ncols,
                                        NULL,
                                        NULL, CUV_ncols,
                                        NULL, CUV_ncols, &work_query, lwork );
            lwork = (int)work_query;
            size_t svdsuperb_nelm = lwork;
            d_work += svdsuperb_nelm;
            ws_needed += svdsuperb_nelm;

            info = LAPACKE_dgesvd_work(  LAPACK_COL_MAJOR,
                                        'A', 'A',
                                        CUV_ncols, CUV_ncols,
                                        _T, ld_T,
                                        d_sigma,
                                        TU, ld_TU,
                                        TV, ld_TV,
                                        svdsuperb,
                                        svdsuperb_nelm);

            double *h_sigma = d_sigma;

            if(rank != 0) { /// truncate according to rank
                finalrank = rank;
                if(rank > size_sigma)
                    finalrank = size_sigma;
            }
            else{ /// truncate according to acc
                int newrank = size_sigma;
                int i;
                for(i=2;i<size_sigma;i++){
                    // ECHO_f(h_sigma[i] )
                    if(h_sigma[i] < relacc)
                    {
                        newrank=i;
                        break;
                    }
                }
                finalrank = newrank; ECHO_I(finalrank)
            }

            //since we store SV we need to scale V by S
            int k;
            for(k = 0; k < finalrank; k++){
                double diagval = h_sigma[k];

                cblas_dscal(CUV_ncols, CBLAS_SADDR(diagval), &TV[k], ld_TV);
            }
        #endif
        Crk[0] = (double)finalrank; ECHO_f(Crk[0])

    //=======================================================================================================
    // construct final U
        #if defined HCORE_GEMM_USE_ORGQR
            double* newUV = d_work;
            size_t newUV_nelm  = _M * finalrank;
            d_work += newUV_nelm;

            info = LAPACKE_dorgqr(  LAPACK_COL_MAJOR,
                                    _M, CUV_ncols, CUV_ncols,
                                    _CU, ld_CU,
                                    qrtauA);
            cblas_dgemm(CblasColMajor,
                        CblasNoTrans, CblasNoTrans,
                        _M, finalrank, CUV_ncols,
                        CBLAS_SADDR(d_one), _CU, ld_CU,
                                            TU, ld_TU,
                        CBLAS_SADDR(d_zero), use_CUV_clone ? _CU_save : newUV, use_CUV_clone ? ld_CU_save : ld_CU);

            if(!use_CUV_clone)
                LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', _M, finalrank, newUV, ld_CU, _CU_save, ld_CU_save);
        #else
            char uplo = 'A';
            int nrows = _M - CUV_ncols;
            int ncols = finalrank;
            LAPACK_dlaset( &uplo, &nrows, &ncols, &d_zero, &d_zero, &(TU[CUV_ncols]), &ld_TU );

            info = LAPACKE_dormqr(  LAPACK_COL_MAJOR,
                                    'L', 'N',
                                    _M, finalrank, CUV_ncols,
                                    _CU, ld_CU,
                                    qrtauA,
                                    TU, ld_TU);

            LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', _M, finalrank, TU, ld_TU, _CU_save, ld_CU_save);
        #endif

    //=======================================================================================================
    // construct final V
        #ifdef HCORE_GEMM_USE_ORGQR
            info = LAPACKE_dorgqr(  LAPACK_COL_MAJOR,
                                    _M, CUV_ncols, CUV_ncols,
                                    _CV, ld_CV,
                                    qrtauB);

            cblas_dgemm(CblasColMajor,
                        CblasNoTrans,
                        #ifdef HCORE_GEMM_USE_KBLAS_ACA
                            CblasNoTrans,
                        #else
                            CblasTrans,
                        #endif
                        _M, finalrank, CUV_ncols,
                        CBLAS_SADDR(d_one), _CV, ld_CV,
                                            TV, ld_TV,
                        CBLAS_SADDR(d_zero), use_CUV_clone ? _CV_save : newUV, use_CUV_clone ? ld_CV_save : ld_CV);

            if(!use_CUV_clone)
                LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', _M, finalrank, newUV, ld_CV, _CV_save, ld_CV_save);
        #else
            #ifdef HCORE_GEMM_USE_KBLAS_ACA
                int TV_pad = CUV_ncols;
                nrows = _M - CUV_ncols;
                ncols = finalrank;
            #else
                int TV_pad = CUV_ncols * ld_TV;
                nrows = finalrank;
                ncols = _M - CUV_ncols;
            #endif
            LAPACK_dlaset( &uplo, &nrows, &ncols, &d_zero, &d_zero, &(TV[TV_pad]), &ld_TV );

            info = LAPACKE_dormqr(  LAPACK_COL_MAJOR,
                                    #ifdef HCORE_GEMM_USE_KBLAS_ACA
                                    'L', 'N',
                                    _M, finalrank, CUV_ncols,
                                    #else
                                    'R', 'T',
                                    finalrank, _M, CUV_ncols,
                                    #endif
                                    _CV, ld_CV,
                                    qrtauB,
                                    TV, ld_TV);

            #ifdef HCORE_GEMM_USE_KBLAS_ACA
                LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'A', _M, finalrank, TV, ld_TV, _CV_save, ld_CV_save);
            #else
                LAPACKE_dge_trans(LAPACK_COL_MAJOR, finalrank, _M, TV, ld_TV, _CV_save, ld_CV_save);
            #endif
        #endif
        ECHO_I(ws_needed);
    ECHO_LN
    }while(0);
}

/***************************************************************************/
// For debugging precision conversion script
    //dormqr
    //zormqr
    //ormqr
    //LAPACKE_dormqr
    //LAPACKE_zormqr    //
    //LAPACKE_ormqr     //
    //dunmqr
    //zunmqr
    //unmqr
    //LAPACKE_dunmqr
    //LAPACKE_zunmqr    //
    //LAPACKE_unmqr     //
