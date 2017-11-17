/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
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
#include <morse.h>
#include <coreblas/cblas.h>
#include <coreblas/lapacke.h>
#include <coreblas/coreblas.h>
#include "timing_zauxiliary.h"


#undef  CBLAS_SADDR                 //FIXME
#define CBLAS_SADDR(_val) (_val)    //FIXME  I should not include this definition
/*--------------------------------------------------------------
 * Check the gemm
 */
double hicma_z_check_gemm(
                   MORSE_enum transA, MORSE_enum transB, int M, int N, int K, //FIXME use z cblas calls for precision generation
                   double alpha, double *A, int LDA,
                   double *B, int LDB,
                   double beta, double *Cmorse,
                   double *Cref, int LDC,
                   double *Cinitnorm, double *Cmorsenorm, double *Clapacknorm )
{
    double beta_const = -1.0;
    double Rnorm;
    double *work = (double *)malloc(chameleon_max(K,chameleon_max(M, N))* sizeof(double));

    *Cinitnorm   = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref,    LDC, work);
    *Cmorsenorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cmorse, LDC, work);

    cblas_dgemm(CblasColMajor, (CBLAS_TRANSPOSE)transA, (CBLAS_TRANSPOSE)transB, M, N, K, //TODO
                CBLAS_SADDR(alpha), A, LDA, B, LDB, CBLAS_SADDR(beta), Cref, LDC);        //TODO

    *Clapacknorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    cblas_daxpy(LDC * N, CBLAS_SADDR(beta_const), Cmorse, 1, Cref, 1);

    Rnorm = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, Cref, LDC, work);

    free(work);

    return Rnorm;
}


