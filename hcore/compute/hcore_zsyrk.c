/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zsyrk.c
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


int syrk_print_index = 0;
extern int use_scratch;

/*
 * N != LDA
 * Leading dim of AU is LDA
 *                                                            M LDA
 *                                                          ______
 *                                                         |  AU  |
 *                                                     K   |______|
 *                                                    ____    
 *                                                   |    |            
 *                                                 N | AV |
 *                                               LDA |    |
 *                                                   |____|            
 *                                       N LDA    X
 *                                      _________
 *                                     |   AV    |
 *         M                        K  |_________|
 *       _____      _____         ____    
 *      |     |    |     |     M |    |               
 *    M | CD  | =  | CD  | + LDA | AU |             
 * LDCD |_____|    |_____|       |____|             
 */
void HCORE_zsyrk(MORSE_enum uplo, MORSE_enum trans,
                int M, int K,
                double alpha,
                const double *AU, int LDAU,
                const double *AV, int LDAV,
                double beta, 
                double *CD, int LDCD,
                double* work
                )
{
    int64_t N = LDAU; //ASSUMPTION FIXME
    int64_t LDA = LDAU;

    /*cblas_zsyrk(*/
            /*CblasColMajor,*/
            /*(CBLAS_UPLO)uplo, (CBLAS_TRANSPOSE)trans,*/
            /*N, N,*/
            /*CBLAS_SADDR(alpha), A, LDA,*/
            /*CBLAS_SADDR(beta), C, LDC);*/

    int64_t A_colfactor_ncols = K;
    /// C = C + alpha * A * A'
    /// C = C + alpha * ( (A^u * (A^v * A^v^T) ) * A^u^T)
    /// A^v * B^v^T
    int64_t bufmtx_nrows = A_colfactor_ncols;
    int64_t bufmtx_ncols = A_colfactor_ncols;
    size_t bufmtx_nelm = bufmtx_nrows * bufmtx_ncols;
    /*size_t bufmtx_nelm = K * K;*/
    //ALLOCATE1 bufmtx_syrk_kk[myid];
    double* bufmtx = NULL;
    if(use_scratch){
        bufmtx = work;
    } else {
        bufmtx = malloc(bufmtx_nelm * sizeof(double));
    }
    assert(bufmtx != NULL);
    double alpha2 = 1.0f;
    double beta2 = 0.0f;
    int64_t ld_AV = M; //ASSUMPTION
    int64_t ld_bufmtx = bufmtx_nrows; //ASSUMPTION

    int64_t AV_nrows = M; //ASSUMPTION
    int64_t AV_ncols = K;
    if(syrk_print_index){
        printf(" SYRK 1 |%d\t|Trans NoTrans AV_ncols:%d AV_ncols:%d AV_nrows:%d alpha2:%.2e ld_AV:%d beta2:%.2e ld_bufmtx:%d AV:%p bufmtx:%p\n", __LINE__,
                AV_ncols, AV_ncols, AV_nrows, alpha2, ld_AV, beta2, ld_bufmtx, AV, bufmtx    
              );
        printf(" SYRK 1 |%d\t|Trans NoTrans K:%d K:%d N:%d alpha2:%.2e LDA:%d beta2:%.2e K:%d AV:%p bufmtx:%p\n", __LINE__,
                                              K, K, N,     alpha2,     LDA,   beta2,     K,   AV, bufmtx
              );
    } 
    /*assert(AV_ncols == K);*/
    /*assert(AV_nrows == N);*/
    /*assert(ld_AV    == LDA);*/
    /*assert(ld_bufmtx== K);*/
    cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, K, K, N, alpha2, AV, LDA, AV, LDA, beta2, bufmtx, K);
    /*bufmtx_nrows = bufmtx_ncols = AV_ncols;  */
    double* AcolBcolT = bufmtx;
    int64_t AcolBcolT_ncols = bufmtx_ncols;

    /// A^u * (A^v * B^v^T)
    int64_t A_rowfactor_nrows = M;
    int64_t A_rowfactor_ncols = K;
    int64_t AU_nrows = M;
    int64_t AU_ncols = K;
    int64_t bufmtx2_nrows = A_rowfactor_nrows; 
    int64_t bufmtx2_ncols = A_rowfactor_ncols;
    //ALLOCATE2
    //bufmtx_syrk_U[myid];
    double* bufmtx2 = NULL;
    
    if(use_scratch){
        bufmtx2 = work + bufmtx_nelm;
    } else {
        /*malloc(bufmtx2_nrows * bufmtx2_ncols * sizeof(double));*/
        bufmtx2 = malloc(M * K * sizeof(double));
    }
    assert(bufmtx2 != NULL);

    int64_t ld_AU = M; //ASSUMPTION
    int64_t ld_bufmtx2 = M; //ASSUMPTION
    if(syrk_print_index){
        printf(" SYRK 2 |%d\t|NoTrans NoTrans AU_nrows:%d AcolBcolT_ncols:%d AU_ncols:%d alpha2:%.2e ld_AU:%d ld_bufmtx:%d beta2:%.2e ld_bufmtx2:%d AU:%p AcolBcolT:%p bufmtx2:%p\n", __LINE__,
                 AU_nrows, AcolBcolT_ncols, AU_ncols, alpha2, ld_AU, ld_bufmtx, beta2, ld_bufmtx2, AU, AcolBcolT, bufmtx2
              );
        printf(" SYRK 2 |%d\t|NoTrans NoTrans M:%d K:%d K:%d alpha2:%.2e LDA:%d K:%d beta2:%.2e M:%d AU:%p AcolBcolT:%p bufmtx2:%p\n", __LINE__,
                                                M,   K,   K,        alpha2, LDA, K,   beta2,     M , AU, AcolBcolT, bufmtx2
              );
        
    } 
    /*assert(AU_nrows           == M);*/
    /*assert(AcolBcolT_ncols    == K);*/
    /*assert(AU_ncols           == K);*/
    /*assert(ld_AU              == LDA);*/
    /*assert(ld_bufmtx2         == M);*/
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, K, K, alpha2, AU, LDA, AcolBcolT, K, beta2, bufmtx2, M);
    bufmtx2_nrows = AU_nrows;
    bufmtx2_ncols = AcolBcolT_ncols;
  

    double* Arow_AcolBcolT = bufmtx2; 
    int64_t ld_C = M; //ASSUMPTION
    /// (A^u * (A^v * B^v^T) ) * B^u^T   
    int64_t Arow_AcolBcolT_nrows = bufmtx2_nrows;
    int64_t Arow_AcolBcolT_ncols = bufmtx2_ncols;
    if(syrk_print_index){
        printf(" SYRK 3 |%d\t|NoTrans Trans Arow_AcolBcolT_nrows:%d AU_nrows:%d Arow_AcolBcolT_ncols:%d alpha:%.2e ld_bufmtx2:%d ld_AU:%d alpha2:%.2e LDCD:%d Arow_AcolBcolT:%p AU:%p CD:%p\n", 
                __LINE__,Arow_AcolBcolT_nrows,  AU_nrows,  Arow_AcolBcolT_ncols,  alpha,  ld_bufmtx2,  ld_AU,  alpha2,  LDCD,  Arow_AcolBcolT,  AU,  CD 
              );
        printf(" SYRK 3 |%d\t|NoTrans Trans M:%d M:%d K:%d alpha:%.2e M:%d LDA:%d alpha2:%.2e LDCD:%d Arow_AcolBcolT:%p AU:%p CD:%p\n", 
                __LINE__, M, M,  K,  alpha,  M,  LDA,  alpha2,  LDCD,  Arow_AcolBcolT,  AU,  CD 
              );
    } 
    /*assert(Arow_AcolBcolT_nrows == M);*/
    /*assert(Arow_AcolBcolT_ncols == K);*/
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 
            M, M, K, 
            alpha, 
            Arow_AcolBcolT, M, 
            AU, LDA, 
            alpha2, 
            CD, LDCD);
    int64_t C_ncols = AU_nrows;
    {
        int i,j;
        for(j = 1; j < M; j++){
            for(i = 0; i < j; i++){
                CD[j*LDCD+i] = beta2;
            }
        }
    }

    if(use_scratch == 0){
        //FREE1
        free(bufmtx);
        //FREE2
        free(bufmtx2);
    }
}


