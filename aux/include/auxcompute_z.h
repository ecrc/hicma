/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file auxcompute_z.h
 *
 * This file contains the declarations of computational auxiliary functions.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#ifndef __AUXCOMPUTE_Z__
#define __AUXCOMPUTE_Z__

#include "hicma_struct.h"
#ifdef MKL
  #include <mkl.h>
  //#pragma message("MKL is used")
#else
  #include <cblas.h>
  #ifdef LAPACKE_UTILS
    #include <lapacke_utils.h>
  #endif
  #include <lapacke.h>
  //#pragma message("MKL is NOT used")
#endif

#include <stdio.h>
#include "morse.h"
#ifndef min
#define min(a, b) ((a) < (b) ? (a) : (b))
#endif


#include "starsh.h"
int HICMA_zuncompress(
        MORSE_enum uplo, MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark);
int HICMA_zuncompress_custom_size(MORSE_enum uplo,
        MORSE_desc_t *AUV, MORSE_desc_t *AD, MORSE_desc_t *Ark,
        int numrows_matrix,
        int numcolumns_matrix,
        int numrows_block,
        int numcolumns_block
        );
int HICMA_zdiag_vec2mat(
        MORSE_desc_t *vec, MORSE_desc_t *mat);
void HICMA_znormest( int M, int N, double *A, double *e, double *work);
void HICMA_zgenerate_problem(
        int probtype, //problem type defined in hicma_constants.h 
        char sym,     // symmetricity of problem: 'N' or 'S'
        double decay, // decay of singular values. Will be used in HICMA_STARSH_PROB_RND. Set 0 for now.
        int _M,       // number of rows/columns of matrix
        int _nb,      // number of rows/columns of a single tile
        int _mt,      // number of tiles in row dimension
        int _nt,      // number of tiles in column dimension
        HICMA_problem_t *hicma_problem // pointer to hicma struct (starsh format will be used to pass coordinate info to number generation and compression phase)
        );
#endif
