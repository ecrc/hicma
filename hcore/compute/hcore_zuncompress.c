/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zuncompress.c
 *
 *  HiCMA HCORE routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"
#include <assert.h>

//FIXME PREVIOUS DECLARION OF CBLAS_SADDR ~/hicma-dev/chameleon/build/include/chameleon/coreblas/include/coreblas.h
#undef CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

int gemmfrk_print_index = 0;
int gemmfrk_print_mat = 0;
int64_t hcfrk_nelm_limit = 289;
void hcfrk_printmat(double * A, int64_t m, int64_t n, int64_t ld){
    printf("M:%d N:%d LD:%d\n[", m, n, ld);
    int64_t i, j, nelm = 0;
    for(i=0;i<m;i++){
        printf("[");
        for(j=0;j<n;j++){
            printf("%+.2e", A[j*ld+i]);
        if(j!=n-1){
            printf(", ");
        }
            //printf("%g ", A[j*tld(descZ)+i]);
            //printf("%g\t", A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
            nelm++;
            if(nelm >= hcfrk_nelm_limit){
                printf("\n");
                return;
            }
        }
        printf("]");
        if(i!=m-1){
            printf(",");
            printf("\n");
        }
    }
    printf("]\n");
}

/*
 * CD=AU*BV'. 
 * Rank of tile AU is in Ark.
 * Rank of tile BV is in Brk.
 * Multiplied tiles must have same rank.
 * CD is dense output.
*/
void HCORE_zuncompress(MORSE_enum transA, int transB,
        int M, int N,
        double alpha, 
        double *AU, 
        double *Ark, 
        int LDA,
        double *BV, 
        double *Brk, 
        int LDB,
        double beta, 
        double *CD, 
        int LDC
)
{
    /*printf("%d|M:%d N:%d K:%d LDA:%d LDB:%d LDC:%d rk:%d acc:%e a:%e b:%e\n",*/
            /*__LINE__, M, N, K, LDA, LDB, LDC, rk, acc, alpha, beta);*/

    /*
     * NOTES:
     * assumptions on matrix dimensions are marked as //ASSUMPTION
     * I am currently allocating and freeing temporary buffers.
     * They are marked as //ALLOCATE
     */
    /*printf("%d %d|%g->%d. %g %g, %g %g\n",  */
            /*__LINE__, __COUNTER__,*/
            /*Crk[0], new_Crk, CU[0], CU[1], CV[0], CV[1]);*/
    //hcore_dgemm(Aij, ik, Ajk, -1, rank, acc); 
    int64_t _Ark = (int64_t)(Ark[0]);
    int64_t _Brk = (int64_t)(Brk[0]);

    if(gemmfrk_print_index){
        printf("Ark:%d Brk:%d M:%d N:%d K:%d  LDA:%d LDB:%d LDC:%d\n",
                _Ark, _Brk, M, N, _Ark, LDA, LDB, LDC  );
    }
    assert(_Ark == _Brk); // ranks must be same
    double one = 1.0, zero = 0.0, minusone = -1.0;
    if(gemmfrk_print_mat){
        hcfrk_printmat(AU,  M, _Ark, LDA); 
        hcfrk_printmat(BV,  N, _Brk, LDB); 
    }
    cblas_dgemm(
            CblasColMajor,
            CblasNoTrans, CblasTrans,
            M, N, _Ark,
            CBLAS_SADDR(one),  AU,  LDA,
            BV,  LDB,
            CBLAS_SADDR(zero), CD, LDC);
    if(gemmfrk_print_mat){
        hcfrk_printmat(CD,  M, N, LDC); 
    }

}






#ifdef ___X
void old(){
    double one = 1.0, zero = 0.0, minusone = -1.0;
    double* CUV = (double*)
        malloc(LDC * LDC * sizeof(double));
//#define DENSEUV
#ifdef DENSEUV
        int64_t inc = 1;
        double* buff = (double*)
            malloc(LDC * LDC * sizeof(double));
        double* buff2 = (double*)
            malloc(LDC * LDC * sizeof(double));
        /*printf("%d AV x BU:%d %d %d\n", */
                /*__LINE__, _Ark, _Brk, K);*/
        cblas_zgemm(
                CblasColMajor,
                CblasTrans, CblasNoTrans,
                _Ark, _Brk, K,
                CBLAS_SADDR(one),  AV,   LDC,
                                   BU,   LDC,
                CBLAS_SADDR(zero), buff, LDC);
        /*printf("%d (AVxBU) x BV %d %d %d\n",*/
                /*__LINE__, _Ark, N, _Brk);*/
        cblas_zgemm(
                CblasColMajor,
                CblasNoTrans, CblasTrans,
                _Ark, N,  _Brk, //ncolsAV can be different from N
                CBLAS_SADDR(one),  buff, LDC,
                                   BV, LDC,
                CBLAS_SADDR(zero), buff2, LDC);
        /*printf("%d AU x ((AVxBU)xBV) %d %d %d\n", */
                /*__LINE__, M, N, _Ark);*/
        if(1){
            cblas_zcopy(LDC*LDC, CD, inc,CUV,inc);
        } else {
            cblas_zgemm(
                    CblasColMajor,
                    CblasNoTrans, CblasTrans,
                    M, N, _Crk,
                    CBLAS_SADDR(one),  CU,  LDC,
                                       CV,  LDC,
                    CBLAS_SADDR(zero), CUV, LDC);
        }
        cblas_zgemm(
                CblasColMajor,
                CblasNoTrans, CblasNoTrans,
                M, N, _Ark,
                CBLAS_SADDR(alpha), AU,    LDC,
                                    buff2, LDC,
                CBLAS_SADDR(beta),  CUV,   LDC);
        /*printf("%d\n", __LINE__);*/

        free(buff);
        free(buff2);
#else
    hcore_dgemm_one(M, N, K, 
            CU, LDC, CV, LDC, _Crk, &new_Crk, 
            AU, LDA, AV, LDA, _Ark, 
            BU, LDB, BV, LDB, _Brk, 
            alpha, rk, acc); 
#endif
    int64_t old_Crk = Crk[0];
    Crk[0] = new_Crk;
    if(0)//AD*BD+=CD
    {
        cblas_zgemm(
                CblasColMajor,
                (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB,
                M, N, K,
                CBLAS_SADDR(alpha), AD, LDA,
                                    BD, LDB,
                CBLAS_SADDR(beta),  CD, LDC);
    } 
    if(0){
        cblas_zgemm(
                CblasColMajor,
                CblasNoTrans, CblasTrans,
                M, N, new_Crk,
                CBLAS_SADDR(one),  CU, LDC,
                                   CV, LDC,
                CBLAS_SADDR(zero), CD, LDC);
    }
//#ifdef aubucheck

    {
        //printf("alpha:%g beta:%g\n",alpha,beta);
#ifndef DENSEUV
        cblas_zgemm(
                CblasColMajor,
                CblasNoTrans, CblasTrans,
                M, N, new_Crk,
                CBLAS_SADDR(one), CU, LDC,
                CV, LDC,
                CBLAS_SADDR(zero), CUV, LDC);
#endif
        int i;
        double dzero = 0.0;
        double dthresh = 1.e-14;
        int diff_nums = 0;
        for(i=0;i<LDC*LDC;i++){
            double diff = fabs((CUV[i] - CD[i])/CUV[i]);
            if(diff > dthresh){
                /*printf("%d: diff:%.2e CD:%.2e CUV:%.2e\n", i,diff, CD[i], CUV[i]);*/
                diff_nums++;
            }
            CD[i] = CUV[i];
        }
        free(CUV);
    /*printf("Ark:%d Brk:%d RANK CHANGE: %d->%d. %g %g, %g %g ndiffs:%d\n",  */
            /*old_Crk, new_Crk, CU[0], CU[1], CV[0], CV[1], diff_nums);*/
    printf("Ark:%d Brk:%d RANK CHANGE: %d->%d\n",  
            _Ark, _Brk, old_Crk, new_Crk);
    /*exit(0);*/
    }
//#endif
}
#endif

/***************************************************************************/

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
