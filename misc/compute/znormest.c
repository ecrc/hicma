/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 * @file znormest.c
 *
 * This file contains the function calculating estimated norm of a matrix. 
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/
#include "morse.h"
#include "auxcompute_z.h"
#include <math.h>

void HICMA_znormest( int M, int N, MORSE_Complex64_t *A, double *e, MORSE_Complex64_t *work)
{
    
    //mkl_dimatcopy('R', 'T', M, N, 1.0, A, N, M);
    
    int i, cnt, maxiter;
    
    MORSE_Complex64_t  e0, alpha, beta;
    MORSE_Complex64_t *X, *SX;
     double tol, normsx, normx;
    X  = work;
    SX = X + M;
    
    /**
     * Let's compute the x vector such that
     *
     *   x_j = sum( A_ij )
     */
    zlange_("I", &M, &N, A, &M, X);
    //LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', M, N, A, M);
    
    /**
     * WARNING: We directly pass the pointer to the full storage of X, since
     * this task depends on "e"
     */
    *e = cblas_dznrm2(N, X, 1);
    
    if ((*e) == 0.0) {
        return;
    } else {
        normx = *e;
    }
    
    alpha = 1.;
    tol = 1.e-1;
    cnt = 0;
    e0  = 0.;
    maxiter = min( 100, N );
    
    while ( (cnt < maxiter) &&
           (fabs((*e) - e0) > (tol * (*e))) )
    {
        e0 = *e;
        
        /**
         * Scale x = x / ||A||
         */
        
        alpha = 1.0/e[0]; i = 1;
        //dscal_(&M, &alpha, X, &i);
        cblas_zscal(M, &alpha, X, i);
        
        /**
         *  Compute Sx = S * x
         */
        e0 = *e; alpha = 1.0; beta = 0.0; i = 1;
        //dgemv_("n", &M, &N, &alpha, A, &M, X, &i, &beta, SX, &i);
        cblas_zgemv(CblasColMajor,CblasNoTrans, M, N, &alpha, A, M, X, i, &beta, SX, i);
        
        /**
         *  Compute x = S' * S * x = S' * Sx
         */
        alpha = 1.0; beta = 0.0; i = 1;
        //dgemv_("t", &M, &N, &alpha, A, &M, SX, &i, &beta, X, &i);
        cblas_zgemv(CblasColMajor,CblasConjTrans, M, N, &alpha, A, M, SX, i, &beta, X, i);
        
        /**
         * Compute ||x||, ||Sx||
         */
        normx  = cblas_dznrm2(M, X,  1);
        normsx = cblas_dznrm2(M, SX, 1);
        
        *e = normx / normsx;
        cnt++;
    }
    
    if ( (cnt >= maxiter) &&
        (fabs((*e) - e0) > (tol * (*e))) ) {
        fprintf(stderr, "mkl_dtwonm_Tile_Async: didn't converge\n");
    }
    
    return;
}

