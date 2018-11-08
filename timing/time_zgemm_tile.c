/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file time_zgemm_tile.c
 *
 * This file shows how to generate tile low-rank (TLR) matrices and perform matrix multiplication.
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
 * The meaning of the descriptors
 * XU : U
 * XV : V
 * XD : U*V
 * X  : the original, non-approximated problem
 * Xrk: rank of U and V, each tile of the matrix is a single number.
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

#include "starpu.h"


#include <assert.h>
#include "hicma_z.h"
#include "auxcompute_z.h"
#include "auxdescutil.h"

#include "hicma.h"
#include "timing_zauxiliary.h"

#include <math.h>
#include <time.h>

#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

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
        fprintf(stdout, "%d:%s\t%d\t%s\t%s\n", myrank, datebuf, __LINE__, __func__, str);\
        fflush(stdout);\
    }
#undef PROGRESS
#define PROGRESS(str)

int store_only_diagonal_tiles = 0;
int global_check = 0;
int global_always_fixed_rank = 0;
int global_fixed_rank = 0;
int num_mpi_ranks;
int diag_nrows = 0;
int main_print_index = 0;
int print_index = 0;
int main_print_mat = 0;
int print_mat = 0;
int use_scratch = 1; // Use scratch memory provided by starpu

double timediff(struct timeval begin, struct timeval end){
    double elapsed = (end.tv_sec - begin.tv_sec) +
        ((end.tv_usec - begin.tv_usec)/1000000.0);
    return elapsed;
}
// This function isused to compare to descriptors in terms of size and numerical values.
// Do not use this function with MPI
// FIXME: This function will be moved to aux/ folder
void check_same(MORSE_desc_t *descL, MORSE_desc_t *descR, char diag, char lower_triangle){
    double *MATL = descL->mat;
    double *MATR = descR->mat;
    if(descL->mb != descR->mb){
        printf("mb: %d %d\n", descL->mb, descR->mb);
    }
    assert(descL->mb == descR->mb);
    if(descL->nb != descR->nb){
        printf("nb: %d %d\n", descL->nb, descR->nb);
    }
    assert(descL->nb == descR->nb);
    if(descL->mt != descR->mt){
        printf("mt: %d %d\n", descL->mt, descR->mt);
    }
    assert(descL->mt == descR->mt);
    if(descL->nt != descR->nt){
        printf("nt: %d %d\n", descL->nt, descR->nt);
    }
    assert(descL->nt == descR->nt);
    int64_t i, j, imt, jnt;
    for(imt=0;imt<descL->mt;imt++){
        for(jnt=0;jnt<descL->nt;jnt++){
            if(diag == 'D' && imt != jnt){
                continue;
            }
            if(lower_triangle == 'L' && imt < jnt){
                continue;
            }
            double *L = &MATL[tsa(descL, imt, jnt)];
            double *R = &MATR[tsa(descR, imt, jnt)];
            for(i=0;i<descL->nb;i++){
                for(j=0;j<descL->nb;j++){
                    double valL = L[j*tld(descL)+i];
                    double valR = R[j*tld(descR)+i];
                    double diff = fabs(valL - valR);
                    double thresh = 1e-14;
                    if(diff > thresh ){
                        printf("Tile:%d,%d. Elm:%d,%d val:%.2e %.2e\n", imt, jnt, i, j, valL, valR);
                        exit(1);
                    }
                    //printf("%g ", A[j*tld(descZ)+i]);
                    //printf("%g\t", A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
                }
            }
        }
    }
}


    int
RunTest(int *iparam, double *dparam, morse_time_t *t_, char* rankfile)
{
    // print progress info only on ROOT process
    if(MORSE_My_Mpi_Rank() != 0)
        print_progress = 0;
    PROGRESS("RunTest started");
    // set alpha and beta which will be used in gemm
    double alpha = 1.0, beta = 1.0;
    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 0;

    HICMA_set_use_fast_hcore_zgemm();

    // Gemm requires more descriptors
    //chameleon/runtime/starpu/control/runtime_descriptor.c
    MORSE_user_tag_size(31,26);
    //MORSE_user_tag_size(31,28);

    global_always_fixed_rank = iparam[IPARAM_HICMA_ALWAYS_FIXED_RANK];

    // get parameters coming from command line
    PASTE_CODE_IPARAM_LOCALS( iparam );
    // set global variable so that p.. files can fill dense matrix
    global_check = check;
    // calculate total number of mpi processes (it is not used for now)
    num_mpi_ranks = P*Q;
    print_index = iparam[IPARAM_HICMA_PRINTINDEX];
    print_mat   = iparam[IPARAM_HICMA_PRINTMAT];
    int64_t _nb = iparam[IPARAM_NB];
    LDB = chameleon_max(K, iparam[IPARAM_LDB]);
    LDC = chameleon_max(M, iparam[IPARAM_LDC]);
    int hicma_maxrank  = iparam[IPARAM_HICMA_MAXRANK];

    // Idealy, in gemm C=AxB,
    // A is M by K
    // B is K by N
    // C is M by N
    // However work with square matrices for now. FIXME
    assert(M==K);
    // Maxrank is N for now FIXME
    assert(M>=N);
    // As a result
    // Dense initial matrices:
    // A is M by M
    // B is M by M
    // C is M by M
    //
    // Low rank matrices:
    // A is M by N
    // B is M by N
    // C is M by N

    int saveN, saveNB;

    saveNB = NB;
    NB = MB;
    /* Allocate Data */
    // initial matrices which will be used for accuracy checking
    // TODO: since HICMA_zgytlr_Tile does not handle NULL dense pointer, have to create them, better condition it with 'check' variable
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, 1, double, MorseRealDouble, LDA, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, 1, double, MorseRealDouble, LDB, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, 1, double, MorseRealDouble, LDC, M, M );
    PROGRESS("desc..'s are allocated");

    // gemm will not have dense diagonals in RND TILED. So remove desc..D's in the future. TODO
    // TODO: since HICMA_zgytlr_Tile does not handle NULL dense pointer, have to create them, better condition it with 'check' variable
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAD, 1, double, MorseRealDouble, LDA, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descBD, 1, double, MorseRealDouble, LDB, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descCD, 1, double, MorseRealDouble, LDC, M, M );
    NB = saveNB;


    //U and V's
    saveN = N;
    N = N * 2;
    saveNB = NB;
    NB = NB * 2;
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUV, 1, double, MorseRealDouble, LDC, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descBUV, 1, double, MorseRealDouble, LDC, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descCUV, 1, double, MorseRealDouble, LDC, M, N );
    N = saveN;
    NB = saveNB;
    PROGRESS("desc..UV's are allocated");

    int mt = descAUV->mt;
    int nt = descAUV->nt;
    assert(mt == nt);

    // rank matrices
    /* tile dimension of rank descriptor must be 1 */
    /* when LD for rk matrices is 1, program exits*/
    int bigMB = MB;
    int bigNB = NB;
    MB = NB = 1;
    int ldrk = mt;//when I use fix number e.g. 2, program exited unexpectedly without any error message
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descArk, 1, double, MorseRealDouble, ldrk, mt, mt);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descBrk, 1, double, MorseRealDouble, ldrk, mt, mt);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descCrk, 1, double, MorseRealDouble, ldrk, mt, mt);
    MB = bigMB;
    NB = bigNB;
    PROGRESS("desc..rk's are allocated");

    int _k = iparam[IPARAM_RK]; //genargs->k
    double _acc = pow(10, -1.0*iparam[IPARAM_ACC]);

    char sym;
    sym = 'N';
    int probtype = iparam[IPARAM_HICMA_STARSH_PROB];
    int maxrank  = iparam[IPARAM_HICMA_STARSH_MAXRANK];
    double ddecay = dparam[IPARAM_HICMA_STARSH_DECAY];
    int  initial_maxrank, final_maxrank;
    double initial_avgrank, final_avgrank;
    HICMA_problem_t hicma_problem;

    //BEGIN: rndtiled
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_RND) {
        hicma_problem.noise    = 0.0;
    }
    //END:   rndtiled

    //BEGIN: geostat
    double theta[3] = {1, 0.1, 0.5};
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_GEOSTAT) {
        hicma_problem.theta = theta;
        hicma_problem.noise = 0.0;
    }
    //END: geostat

    //BEGIN: ss
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_SS) {
        // Correlation length
        hicma_problem.beta  = 0.1;
        //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

        // Smoothing parameter for Matern kernel
        hicma_problem.nu    = 0.5;
        // Shift added to diagonal elements
        hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
        hicma_problem.noise = 5.e-4; //works for 640K
    }
    //END: ss
    PROGRESS("generating coordinates started");
    HICMA_zgenerate_problem(probtype, sym, ddecay, M, MB, mt, nt, &hicma_problem);
    PROGRESS("generating coordinates ended");
    mpiF = hicma_problem.starsh_format; // This is assignment will be hidden from user in release

	// enforce compression of diagonal tiles
	int compress_diag = 1;
    PROGRESS("zgytlr starting");
    //desc[A,B,C] are original problem if global check is enabled
    //desc[A,B,C]D are original problem. Temporary buffers should be used in hcore_zgytlr TODO Do not use desc..D in this gemm code
    HICMA_zgytlr_Tile(MorseUpperLower, descAUV, descAD, descArk, 0, maxrank, _acc, compress_diag,  descA);
    HICMA_zgytlr_Tile(MorseUpperLower, descBUV, descBD, descBrk, 0, maxrank, _acc, compress_diag,  descB);
    HICMA_zgytlr_Tile(MorseUpperLower, descCUV, descCD, descCrk, 0, maxrank, _acc, compress_diag,  descC);
    PROGRESS("zgytlr finished");

    int set_diag = 0;
    double diagVal = M;

    /* Save A for check */
    PROGRESS("pasting original dense descC into C2 started");
    // Adense: original dense problem.
    PASTE_TILE_TO_LAPACK( descC, C2, check, double, LDC, M );
    PROGRESS("pasting original dense descC into C2 finished");

    PROGRESS("gemm started");
    START_TIMING();
    HICMA_zgemm_Tile( MorseNoTrans, MorseNoTrans, alpha,  //TODO
            descAUV, descArk,
            descBUV, descBrk,
            beta,
            descCUV, descCrk, _k, maxrank, _acc
            );
    STOP_TIMING();
    fflush(stderr);
    fflush(stdout);
    PROGRESS("gemm finished");

    if(check){

        PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, M );
        PASTE_TILE_TO_LAPACK( descB, B, check, double, LDB, M );
        HICMA_zuncompress(MorseUpperLower, descCUV, descC, descCrk);

        PASTE_TILE_TO_LAPACK( descC, C, check, double, LDC, M );


        dparam[IPARAM_RES] = hicma_z_check_gemm( MorseNoTrans, MorseNoTrans, M, M, M,
                alpha, A, LDA, B, LDB, beta, C, C2, LDC,
                &(dparam[IPARAM_ANORM]),
                &(dparam[IPARAM_BNORM]),
                &(dparam[IPARAM_XNORM]));
        free(A); free(B); free(C); free(C2);
    }

    PASTE_CODE_FREE_MATRIX( descAUV );
    PASTE_CODE_FREE_MATRIX( descBUV );
    PASTE_CODE_FREE_MATRIX( descCUV );
    PROGRESS("desc..UV's are freed");

    PASTE_CODE_FREE_MATRIX( descAD );
    PASTE_CODE_FREE_MATRIX( descBD );
    PASTE_CODE_FREE_MATRIX( descCD );
    PROGRESS("desc..D's are freed");

    PASTE_CODE_FREE_MATRIX( descArk );
    PASTE_CODE_FREE_MATRIX( descBrk );
    PASTE_CODE_FREE_MATRIX( descCrk );
    PROGRESS("desc..rk's are freed");

    PASTE_CODE_FREE_MATRIX( descA );
    PASTE_CODE_FREE_MATRIX( descB );
    PASTE_CODE_FREE_MATRIX( descC );
    PROGRESS("desc..'s are freed");
    PROGRESS("freed descs");
    return 0;
}
