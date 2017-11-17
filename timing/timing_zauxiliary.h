/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file timing_zauxiliary.h
 *
 * This file contains the declarations of auxiliary functions used for timing experiments.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
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

#ifndef TIMING_ZAUXILIARY_H
#define TIMING_ZAUXILIARY_H

double hicma_z_check_gemm(MORSE_enum transA, MORSE_enum transB, int M, int N, int K,
                   double alpha, double *A, int LDA,
                   double *B, int LDB,
                   double beta, double *Cmorse,
                   double *Cref, int LDC,
                   double *Cinitnorm, double *Cmorsenorm, double *Clapacknorm );



#endif /* TIMING_ZAUXILIARY_H */
