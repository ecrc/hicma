/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *           ?          All rights reserved.
 **/
/**
 * @file time_zpotrf_tile.c
 *
 * This file shows how to generate tile low-rank (TLR) matrix and factorize it using Cholesky factorization.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

/**
 * The meaning of the descriptors:
 * - AUV: U and V, side by side
 * - AD : U*V
 * - A  : the original, non-approximated problem
 * - Ark: rank of U and V, each tile of the matrix is a single integer in fact a double.
 *
 **/

#include "morse.h"
#include "timing.h"
#include "hicma_constants.h"
#include "hicma_struct.h"
#include "hicma_z.h"
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
//#include <mpi.h> //MPI_Wtime()

#include "starpu.h"
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

#include "starsh-spatial.h"

#include <assert.h>
#include "hicma_z.h"
#include "auxcompute_z.h"
#include "auxdescutil.h"
#include "hicma.h"
#include <math.h>
#include <time.h>
#include"hicma_common.h"
#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

char norm = 'F';
// zgytlr uses starsh in MPI mode.
STARSH_blrf *mpiF;

int print_progress = 1;   // Print progress about the execution
char datebuf[128];
time_t timer;
struct tm* tm_info;
#define PROGRESS(str) \
    if(print_progress){ \
        int myrank = MORSE_My_Mpi_Rank();\
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S",tm_info); \
        fprintf(stderr, "%d:%s\t%d\t%s\t%s\n", myrank, datebuf, __LINE__, __func__, str);\
        fflush(stderr);\
    }
//#undef PROGRESS
//#define PROGRESS(str)

int store_only_diagonal_tiles = 0;
int global_check = 0;
int global_always_fixed_rank = 0;
int global_fixed_rank = 0;
int global_omit_computation = 1;
int num_mpi_ranks;
int run_getrf = 1;
int diag_nrows = 0;
int main_print_index = 0;
int print_index = 0;
int print_index_end = 0;
int main_print_mat = 0;
int print_mat = 0;
int use_scratch = 1; // Use scratch memory provided by starpu
int calc_rank_stat = 1; 
int isPrecentageused = 1;
float precent1 = 0.75;
float precent2 = 0.5;
double timediff(struct timeval begin, struct timeval end){
    double elapsed = (end.tv_sec - begin.tv_sec) +
        ((end.tv_usec - begin.tv_usec)/1000000.0);
    return elapsed;
}

void Acoustic_Init(int *nip, int *ntrian);
void AL4SAN_RHS( double _Complex *rhs ,int *nip, int *ntrian);

#if defined(HICMA_COMPLEX)
    int
RunTest(int *iparam, double *dparam, morse_time_t *t_, char* rankfile)
{
    PASTE_CODE_IPARAM_LOCALS( iparam );  
    //HICMA_set_print_mat();
    /*HICMA_set_print_index();*/
    printf("\n ************************I am in complex part ************************\n");
    PROGRESS("RunTest started");

    MORSE_user_tag_size(31,27);

    NB = MB;
    N = M;
    int global_nt=N/NB;
    int local_nt=NB/nip;
    printf("\n nip:%d , ntrian:%d N:%d, NB:%d\n", nip, ntrian, N, NB);
    printf("\n global_nt:%d local_nt:%d N:%d NB:%d \n", global_nt, local_nt, N, NB);
    Acoustic_Init(&nip, &ntrian);

    PASTE_CODE_ALLOCATE_MATRIX_TILE( descDense, 1, MORSE_Complex64_t, MorseComplexDouble, (nip* ntrian), (nip* ntrian), (nip* ntrian) ); 
    PROGRESS("descDense is allocated");

    PROGRESS("HICMA_zgene started");
    struct timeval tvalBefore, tvalAfter;
    gettimeofday (&tvalBefore, NULL);
    HICMA_zgene_Tile( ntrian,  nip,  local_nt,  NB, descDense);
    gettimeofday (&tvalAfter, NULL);
    if(MORSE_My_Mpi_Rank()==0){
        printf("Tgenerate:%g\n", 
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
              );
        fflush(stderr);
        fflush(stdout);
    }
    if(0)printf("%s %d: descACO LDA:%d M:%d MB:%d NB:%d\n", __FILE__, __LINE__, (nip* ntrian), (nip* ntrian), MB, NB); fflush(stdout);
    PROGRESS("HICMA_zgene ended");
    num_mpi_ranks = P*Q;
    starpu_task_wait_for_all();     
    START_TIMING();
    MORSE_zgetrf_nopiv_Tile( descDense);
    PASTE_CODE_FREE_MATRIX( descDense );
    STOP_TIMING();
    return 0;
}

#endif
