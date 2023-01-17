/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * This file contains auxilary functions for routines in timing folder.
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
// #include "libhqr.h"
#include <hicma.h>
#include <coreblas/hicma_cblas.h>
#include "coreblas/hicma_lapacke.h"
#include <coreblas/hicma_coreblas.h>
#include "timing_dauxiliary.h"


#undef  CBLAS_SADDR                 //FIXME
#define CBLAS_SADDR(_val) (_val)    //FIXME  I should not include this definition
/*--------------------------------------------------------------
 * Check the gemm
 */
double hicma_d_check_gemm(
                   HICMA_enum transA, HICMA_enum transB, int M, int N, int K, //FIXME use z cblas calls for precision generation
                   double alpha, double *A, int LDA,
                   double *B, int LDB,
                   double beta, double *Chicma,
                   double *Cref, int LDC,
                   double *Cinitnorm, double *Chicmanorm, double *Clapacknorm )
{
    double beta_const = -1.0;
    double Rnorm;
    double *work = (double *)malloc(hicma_max(K,hicma_max(M, N))* sizeof(double));

    *Cinitnorm   = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref,    LDC, work);
    *Chicmanorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Chicma, LDC, work);

    cblas_dgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, //TODO
                CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC);        //TODO

    *Clapacknorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    cblas_daxpy(LDC * N, CBLAS_SADDR(beta_const), Chicma, 1, Cref, 1);

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    free(work);

    return Rnorm;
}


