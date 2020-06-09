/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zgytlr.c
 *
 *  HiCMA HCORE routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
//#include "hcore/include/hcore.h"
#include "morse.h"
#include "hcore_z.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday
#include <stdlib.h>//FIXME for malloc

#define Stringize( L )     #L 
#define MakeString( M, L ) M(L)
#define $Line MakeString( Stringize, __LINE__ )
#define Reminder __FILE__ "(" $Line ") : Reminder: "

#define COMPLEX
#undef REAL
//  HCORE_zgytlr - Generate a tile for random matrix.

#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"
#include "hicma_common.h"
extern STARSH_blrf *mpiF;
extern int print_index;
int print_index_end;
extern int store_only_diagonal_tiles;
extern int global_check;
extern int print_mat;

int use_rsdd = 1; 
int zgytlr_print_index = 0;
int zgytlr_print_index_end = 0;


extern void _printmat_complex(MORSE_Complex64_t * A, int m, int n, int ld);

//extern int global_always_fixed_rank; //FIXME this does not work, I do not know why
//extern int global_fixed_rank;
int global_always_fixed_rank;
int global_fixed_rank;
int gytlr_tile_ii = -1;
int gytlr_tile_jj = -1;


int dense_dsvfr(int size, double *S, double tol)
//! Returns rank of double precision singular values.
/*! Tries ranks `size`, `size`-1, `size`-2 and so on. May be accelerated by
 * binary search, but it requires additional temporary memory to be allocated.
 *
 * @param[in] size: Number of singular values.
 * @param[in] S: Array of singular values.
 * @param[in] tol: Relative error tolerance.
 * @return rank in terms of relative error in Frobenius norm.
 * */
{
    int i;
    double err_tol = 0;
    // Compute Frobenius norm by `S`
    for(i = 0; i < size; i++)
        err_tol += S[i]*S[i];
    // If all elements of S are zeros, then rank is 0
    if(err_tol == 0)
        return 0;
    // If Frobenius norm is not zero, then set rank as maximum possible value
    i = size;
    // Set tolerance
    err_tol *= tol*tol;
    double tmp_norm = S[size-1]*S[size-1];
    // Check each possible rank
    while(i > 1 && err_tol >= tmp_norm)
    {
        i--;
        tmp_norm += S[i-1]*S[i-1];
    }
    return i;
}
double* dense_dlrrsdd(int nrows, int ncols, MORSE_Complex64_t  *D, int ldD, MORSE_Complex64_t  *U,
        int ldU, MORSE_Complex64_t  *V, int ldV, int *rank, int maxrank, int oversample,
        double tol, MORSE_Complex64_t  *work, int lwork, int *iwork)
//! Randomized SVD approximation of a dense double precision matrix.
/*! This function calls LAPACK and BLAS routines, so integer types are int
 * instead of @ref STARSH_int.
 *
 * @param[in] nrows: Number of rows of a matrix.
 * @param[in] ncols: Number of columns of a matrix.
 * @param[in,out] D: Pointer to dense matrix.
 * @param[in] ldD: leading dimensions of `D`.
 * @param[out] U: Pointer to low-rank factor `U`.
 * @param[in] ldU: leading dimensions of `U`.
 * @param[out] V: Pointer to low-rank factor `V`.
 * @param[in] ldV: leading dimensions of `V`.
 * @param[out] rank: Address of rank variable.
 * @param[in] maxrank: Maximum possible rank.
 * @param[in] oversample: Size of oversampling subset.
 * @param[in] tol: Relative error for approximation.
 * @param[in] work: Working array.
 * @param[in] lwork: Size of `work` array.
 * @param[in] iwork: Temporary integer array.
 * */
{
    int mn = nrows < ncols ? nrows : ncols;
    int mn2 = maxrank+oversample;
    int i;
    if(mn2 > mn)
        mn2 = mn;
    //size_t svdqr_lwork = (4*mn2+7)*mn2;
    //if(svdqr_lwork < ncols)
    //    svdqr_lwork = ncols;
    MORSE_Complex64_t *X, *Q, *tau, *svd_U, *svd_V, *svdqr_work; 
    double *svd_S, *svd_rwork;

    int mx = nrows>ncols?nrows:ncols;
    size_t s1 = 5*mn*mn + 5*mn;
    size_t s2 = 2*mx*mn + 2*mn*mn + mn;
    size_t lrwork = s1>s2?s1:s2;
    svd_rwork = malloc(sizeof(*svd_rwork) * lrwork);
    //double *X, *Q, *tau, *svd_U, *svd_S, *svd_V, *svdqr_work;
    X = work;
    Q = X+(size_t)ncols*mn2;
    svd_U = Q+(size_t)nrows*mn2;
    svd_S = svd_U+(size_t)mn2*mn2;
    tau = svd_S;
    svd_V = svd_S+mn2;
    svdqr_work = svd_V+ncols*mn2;
    int svdqr_lwork = lwork-(size_t)mn2*(2*ncols+nrows+mn2+1);
    int iseed[4] = {0, 0, 0, 1};
    MORSE_Complex64_t zero = (MORSE_Complex64_t) 0.0;
    MORSE_Complex64_t one = (MORSE_Complex64_t) 1.0;
    // Generate random matrix X
    LAPACKE_zlarnv_work(3, iseed, nrows*mn2, X);
    // Multiply by random matrix
    cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows, mn2,
            ncols, &one, D, ldD, X, ncols, &zero, Q, nrows);
    // Get Q factor of QR factorization
    LAPACKE_zgeqrf_work(LAPACK_COL_MAJOR, nrows, mn2, Q, nrows, tau,
            svdqr_work, svdqr_lwork);
    LAPACKE_zungqr_work(LAPACK_COL_MAJOR, nrows, mn2, mn2, Q, nrows, tau,
            svdqr_work, svdqr_lwork);
    // Multiply Q by initial matrix
    cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, mn2, ncols,
            nrows, &one, Q, nrows, D, ldD, &zero, X, mn2);
    // Get SVD of result to reduce rank
    int info = LAPACKE_zgesdd_work(LAPACK_COL_MAJOR, 'S', mn2, ncols, X, mn2,
            svd_S, svd_U, mn2, svd_V, mn2, svdqr_work, svdqr_lwork, svd_rwork, iwork);
    if(info != 0) {
        printf("%s %d: LAPACKE_dgesdd_work info=%d", __FILE__, __LINE__, info);
        exit(-1);
    }
    // Get rank, corresponding to given error tolerance
    *rank = dense_dsvfr(mn2, svd_S, tol);
    //printf("%s %d: info:%d rank:%d\n", __FILE__, __LINE__, info, *rank);
    if(info == 0 && *rank <= maxrank)
    // If far-field block is low-rank
    {
        cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, nrows, *rank,
                mn2, &one, Q, nrows, svd_U, mn2, &zero, U, ldU);
        for(i = 0; i < *rank; i++)
        {
            cblas_zcopy(ncols, svd_V+i, mn2, V+i*(size_t)ldV, 1);
            MORSE_Complex64_t __svd_Si = (MORSE_Complex64_t) svd_S[i];
            cblas_zscal(ncols, &__svd_Si, V+i*(size_t)ldV, 1);
        }
    }
    else
    // If far-field block is dense, although it was initially assumed
    // to be low-rank. Let denote such a block as false far-field block
        *rank = -1;
    free(svd_rwork);
    return svd_S;
}

void HCORE_zgytlr( int m, int n, /*dimension of squareAD*/
        MORSE_Complex64_t *AU,
        MORSE_Complex64_t *AV,
        MORSE_Complex64_t *AD,
        double *Ark,
        int lda,
        int ldu,
        int ldv,
        int bigM, int ii, int jj, unsigned long long int seed,
        int maxrank, double tol, int compress_diag,
        MORSE_Complex64_t *Dense
        )
{
    int64_t i, j;
    //printf("m:%d n:%d bigM:%d m0:%d n0:%d\n", m, n, bigM, m0, n0);
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    if(Ark[0] >= 1 && Ark[0]<maxrank) {
        maxrank = (int)Ark[0];
    }
    if(print_index || zgytlr_print_index){
        fprintf(stderr, "%d+GYTLR\t|(%d,%d) irk:%g mxrk:%d m:%d n:%d lda:%d ldu:%d ldv:%d\n",MORSE_My_Mpi_Rank(), ii, jj, Ark[0], maxrank, m, n, lda, ldu, ldv);
    }

    int shape[2];
    int rank = 0;

    if(isPrecentageused){
        int sp=precent2*rank;
        if (sp==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); sp=1;}
        size_t nelem=(ldu*rank)+(ldu*sp);
        //printf("\n nelem:%d, rank:%d\n", nelem, rank);
        AV=&(AU[nelem]);
    }

    if(print_mat){
        printf("%d\tgytlr-UV-input\n", __LINE__);
        //_printmat_complex(AD, m, n, lda);
        //_printmat_complex(AU, m, rank, ldu);
        //_printmat_complex(AV, ldv, rank, ldv);
        _printmat_complex(Dense, m, n, lda);
    }
    char chall = 'A';
    if(ii == jj){ /** copy from Dense to AD */
        LAPACK_zlacpy(&chall,
                &m, &n, Dense, &lda, AD, &lda);
        return;
    }

    //starsh_dense_zlrrsdd(m, n, AD, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);
    //lapack_int LAPACKE_zgesvd( int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, double* s, lapack_complex_double* u, lapack_int ldu, lapack_complex_double* vt, lapack_int ldvt, double* superb );
    //https://software.intel.com/en-us/node/521150
    struct timeval tvalBefore_svd, tvalAfter_svd;  // removed comma
    gettimeofday (&tvalBefore_svd, NULL);
    int info = 0;
    if(use_rsdd) {
        MORSE_Complex64_t  *work;
        int *iwork;
        int oversample = 10;
        int mn = m;
        int mn2 = maxrank+oversample;
        if(mn2 > mn)
            mn2 = mn;
        // Get size of temporary arrays
        size_t lwork = n, lwork_sdd = (4*mn2+7)*mn2;
        if(lwork_sdd > lwork)
            lwork = lwork_sdd;
        lwork += (size_t)mn2*(2*n+m+mn2+1);
        size_t liwork = 8*mn2;
        // Allocate temporary arrays
        //STARSH_MALLOC(iwork, liwork);
        iwork = malloc(sizeof(*iwork) * liwork);
        if(iwork == NULL) {
            fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! liwork:%d", __FILE__, __func__, __LINE__, liwork);
            exit(-1);
        }
        //STARSH_MALLOC(work, lwork);
        work = malloc(sizeof(*work) * lwork);
        if(work == NULL) {
            fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! lwork:%d", __FILE__, __func__, __LINE__, lwork);
            exit(-1);
        }
        //AD is m x n. AU and AV are m x maxrank and n x maxrank correspondingly
        //starsh_kernel_drsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        //starsh_dense_dlrrsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        double* svd_S = dense_dlrrsdd(m, n, Dense, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);

        if(0){ /** Print singular values */
            char* str = calloc(10*n, sizeof(char));
            char* pstr = str;
            //pstr+=sprintf(pstr, "ii:%d jj:%d m:%d n:%d\n", ii, jj, m, n);
            int limsv=m;
            pstr+=sprintf(pstr, "%d,%d,%d:", ii, jj, limsv);
            for(int i = 0; i < limsv; i++){
                pstr+=sprintf(pstr, "%.2e ", svd_S[i]);
            }
            pstr+=sprintf(pstr, "\n");
            printf("%s", str);
            free(str);
        }

        free(work);
        free(iwork);

    } else {
        MORSE_Complex64_t* Dense_clone = calloc(m * lda, sizeof(MORSE_Complex64_t));
        LAPACK_zlacpy(&chall,
                &m, &n, Dense, &lda, Dense_clone, &lda);


        int nsv = m < n ? m : n;
        double* sigma = calloc(nsv, sizeof(double));
        double* svdsuperb = calloc(nsv, sizeof(double));
        MORSE_Complex64_t *__AV = calloc(nsv * ldv, sizeof(MORSE_Complex64_t));
        info = LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, Dense_clone, lda, sigma, AU, ldu, __AV, ldv, svdsuperb);
        free(Dense_clone);
        //if(info == 0) the execution is successful.
        if(info < 0) { 
            printf("the %d-th parameter had an illegal value\n", info);
        }
        if(info > 0) {
            printf("if ?bdsqr did not converge, %d specifies how many superdiagonals of the intermediate bidiagonal form B did not converge to zero (see the description of the superb parameter for details).\n", info);
        }
        if(info != 0){
            fprintf(stderr,
                    "%s %d ERROR in LAPACKE_zgesvd() info=%d "
                    "1:m=%d, 2:n=%d, 3:Dense=%p, 4:lda=%d, 5:sigma=%p,"
                    "6:U=%p, 7:ldu=%d, 8:V=%p, 9:ldv:%d,"
                    "10:svdsuperb:%p"
                    "\n",
                    __FILE__, __LINE__, info,
                    m, n, Dense, lda, sigma,
                    AU, ldu, __AV, ldv,
                    svdsuperb);
            exit(-1);
        }

        for(i=0; i<nsv; i++){
            if(sigma[i] < tol){
                break;
            }
        }
        rank = i;
        if(print_mat){
            printf("Singular values: ");
            int i;
            for(i=0; i<nsv; i++){
                printf("%.2e ", sigma[i]);
            }
            printf("\n");
        }
        if(rank == nsv){ //means that tile is dense.
            rank = nsv;
            fprintf(stderr, "%s %s %d: Dense off-diagonal block (%d,%d)\n", __FILE__, __func__, __LINE__, ii, jj);
            exit(0);
        }
        int k;
        for(k = 0; k < rank; k++){
            double ddiagval = sigma[k];
            MORSE_Complex64_t zdiagval = ddiagval+0.0 * _Complex_I;
            /*printf("k:%d %+.4f%+.4fi\n", k, creal(zdiagval), cimag(zdiagval));*/
            cblas_zscal(n, CBLAS_SADDR(zdiagval), &__AV[k], ldv);
        } 
        LAPACKE_zge_trans(LAPACK_COL_MAJOR, rank, n,
                __AV, ldv, AV, ldv);

        free(sigma);
        free(svdsuperb);
        free(__AV);
    }
    gettimeofday (&tvalAfter_svd, NULL);
    double time_svd = 
                (tvalAfter_svd.tv_sec - tvalBefore_svd.tv_sec)
                +(tvalAfter_svd.tv_usec - tvalBefore_svd.tv_usec)/1000000.0;
    Ark[0] = rank;
    if(print_mat){
        printf("%d\tgytlr-UV-output\n", __LINE__);
        /*_printmat_complex(AD, m, n, lda);*/
        _printmat_complex(AU, m, rank, ldu);
        _printmat_complex(AV, n, rank, ldv);
        /*_printmat_complex(Dense, m, n, lda);*/
    }
    if(print_index || print_index_end || zgytlr_print_index_end){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-GYTLR\t|(%d,%d) rk:%g m:%d n:%d lda:%d ldu:%d ldv:%d\t\t\t\t\tGYTLR: %.4f svd:%.4f\n",MORSE_My_Mpi_Rank(),ii,jj,Ark[0],m, n, lda, ldu, ldv,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0,
                time_svd
               );
    }
    return;
}




