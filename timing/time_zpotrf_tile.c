/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
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
int run_potrf = 1;
int diag_nrows = 0;
int main_print_index = 0;
int print_index = 0;
int print_index_end = 0;
int main_print_mat = 0;
int print_mat = 0;
int use_scratch = 1; // Use scratch memory provided by starpu
int calc_rank_stat = 1; 

void fwrite_array(int m, int n, int ld, double* arr, char* file){
        FILE* fp = fopen(file, "w");
        if(fp == NULL){
            fprintf(stderr, "File %s cannot be opened to write\n", file);
            exit(1);
        }
        int i, j;
        fprintf(fp, "%d %d\n", m, n);
        for(i = 0; i < m; i++){
            for(j = 0; j < n; j++){
                fprintf(fp, "%d\t", (int)arr[ld*j+i] );
            }
            fprintf(fp, "\n" );
        }
        fclose(fp);
}

double timediff(struct timeval begin, struct timeval end){
    double elapsed = (end.tv_sec - begin.tv_sec) +
        ((end.tv_usec - begin.tv_usec)/1000000.0);
    return elapsed;
}



    int
RunTest(int *iparam, double *dparam, morse_time_t *t_, char* rankfile)
{
    // print progress info only on ROOT process
    if(MORSE_My_Mpi_Rank() != 0)
        print_progress = 0;
    PROGRESS("RunTest started");

    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 1;
    //chameleon/runtime/starpu/control/runtime_descriptor.c
    //MORSE_user_tag_size(31,26);
    //MORSE_user_tag_size(31,29);
    MORSE_user_tag_size(31,27);// When I added tile_to_lapack for descArk, I got not enough number of desc error

    // get parameters coming from command line
    PASTE_CODE_IPARAM_LOCALS( iparam );

    // set global variable so that p.. files can fill dense matrix
    global_check = check;
    // calculate total number of mpi processes (it is not used for now)
    num_mpi_ranks = P*Q;
    print_index = iparam[IPARAM_HICMA_PRINTINDEX];
    print_index_end = iparam[IPARAM_HICMA_PRINTINDEXEND];
    print_mat   = iparam[IPARAM_HICMA_PRINTMAT];
    int64_t _nb = iparam[IPARAM_NB];
    LDA = chameleon_max(M, iparam[IPARAM_LDA]);
    int hicma_maxrank  = iparam[IPARAM_HICMA_MAXRANK];
    global_always_fixed_rank = iparam[IPARAM_HICMA_ALWAYS_FIXED_RANK];

    int saveNB = NB;
    NB = MB;
    size_t ncols_AD;
    int saveP = P;
    int saveQ = Q;
    if (store_only_diagonal_tiles == 1) {
        ncols_AD = MB;
    } else {
        ncols_AD = M;
    }
    int saveN = N;
    N = ncols_AD;
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAD, 1, double, MorseRealDouble, LDA, M, N );
    N = saveN;
    P = saveP;
    Q = saveQ;
    PROGRESS("descAD is allocated");

    size_t ncols_Dense;
    size_t ld_Dense;
    int saveMB = MB;
    if(check == 0) {
        ncols_Dense = MT;
        MB = NB = 1;
        ld_Dense = MT;
    } else {
        ncols_Dense = M;
        ld_Dense = M;
    }
    /*descDense is full matrix if numerical accuracy will be checked.
     * Otherwise it is MB-by-MB matrix with 1-by-1 tiles */
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descDense, 1, double, MorseRealDouble, ld_Dense, ncols_Dense, ncols_Dense );
    if(check == 0) {
        MB = saveMB;
    } else {
    }
    PROGRESS("descDense is allocated");
    NB = saveNB;

    int MTMB = MT * MB; // roundup number of rows/columns for AUV 
    int nrows_AUV = MTMB;
    int ld_AUV = MTMB;
    // allocate descUV
    saveN = N;
    N = N * 2;
    saveNB = NB;
    NB = NB * 2;
    //printf("N:%d NB:%d\n", N, NB);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUV, 1, double, MorseRealDouble, ld_AUV, nrows_AUV, N );
    N = saveN;
    NB = saveNB;
    PROGRESS("descAUV is allocated");

    /* tile dimension of rank descriptor must be 1 */
    /* when LD for rk matrices is 1, program exits*/
    int bigMB = MB;
    int bigNB = NB;
    MB = NB = 1;
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descArk, 1, double, MorseRealDouble, MT, MT, NT);
    PROGRESS("descA's are allocated");
    MB = bigMB;
    NB = bigNB;

    int diag_dense = 1;
    int fixedrank = iparam[IPARAM_RK]; //genargs->k
    double fixedacc = pow(10, -1.0*iparam[IPARAM_ACC]);

    char sym;
    if (run_potrf)
        sym = 'S';
    else
        sym = 'N';
    int probtype = iparam[IPARAM_HICMA_STARSH_PROB];
    int maxrank  = iparam[IPARAM_HICMA_STARSH_MAXRANK];
    //double ddecay = pow(10, -1.0*iparam[IPARAM_HICMA_STARSH_DECAY]);
    double ddecay = dparam[IPARAM_HICMA_STARSH_DECAY];


    int  initial_maxrank, final_maxrank;
    double initial_avgrank, final_avgrank;
    HICMA_problem_t hicma_problem;

    hicma_problem.ndim = 2;

    //BEGIN: rndtiled
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_RND) {
        hicma_problem.noise    = 1.0; //value added to diagonal
    }
    //END:   rndtiled

    //BEGIN: geostat
    //double theta[3] = {1, 0.1, 0.5}; //initially
    double theta[3] = {
        1.0, //sigma 
        0.01, //beta
        10.0 //nu Aleks used 10.0 in his paper
    };
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_GEOSTAT) {
        hicma_problem.theta = theta;
        hicma_problem.noise = 0.0;
        hicma_problem.noise = 1.e-2;
        hicma_problem.kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
    }
    //END: geostat

    //BEGIN: ss
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_SS) {
        //sigma=1.0 default value line 193 of stars-h/src/applications/spatial.c 
        // Correlation length
        hicma_problem.beta  = 0.1;
        //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

        // Smoothing parameter for Matern kernel
        hicma_problem.nu    = 0.5;
        // Shift added to diagonal elements
        hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
        hicma_problem.noise = 5.e-4; //works for 640K but does not work for 10M
        hicma_problem.noise = 1.e-2; //
    }
    //END: ss
    
    //BEGIN: edsin
    if(iparam[IPARAM_HICMA_STARSH_PROB] == HICMA_STARSH_PROB_EDSIN) {
        // Wave number, >= 0
        hicma_problem.wave_k = dparam[IPARAM_HICMA_STARSH_WAVE_K];
        hicma_problem.diag = M; 
    //printf("%s %d: %g\n", __FILE__, __LINE__, hicma_problem.wave_k);
    }
    //END: edsin

    PROGRESS("generating coordinates started");
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    HICMA_zgenerate_problem(probtype, sym, ddecay, M, MB, MT, NT, &hicma_problem);
    gettimeofday (&tvalAfter, NULL);
    if(MORSE_My_Mpi_Rank()==0){
        printf("Tproblem:%g\n",
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
              );
        fflush(stderr);
        fflush(stdout);
    }
    PROGRESS("generating coordinates ended");
    mpiF = hicma_problem.starsh_format; // This is assignment will be hidden from user in release

	// DO NOT enforce compression of diagonal tiles
	int compress_diag = 0;
    PROGRESS("nompi zgytlr starting");
    //descDense original problem
    gettimeofday (&tvalBefore, NULL);
    HICMA_zgytlr_Tile(MorseLower, descAUV, descAD, descArk, 0, maxrank, fixedacc, compress_diag, descDense);
    gettimeofday (&tvalAfter, NULL);
    if(MORSE_My_Mpi_Rank()==0){
        printf("Tcompress:%g\n", 
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
        fflush(stderr);
        fflush(stdout);
    }
    PROGRESS("nompi zgytlr finished");
    fflush(stderr);
    fflush(stdout);
    /*return 0; //TODO*/

    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descArk, Ark_initial, 1, double, MT, NT );
        if(MORSE_My_Mpi_Rank()==0){

            sprintf(rankfile, "%s-1", rankfile);
            fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);

            HICMA_stat_t hicma_statrk_initial;
            zget_stat(MorseLower, Ark_initial, MT, NT, MT,  &hicma_statrk_initial);
            printf("initial_ranks:");
            zprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }

    if (global_always_fixed_rank == 1) {
        fprintf(stderr, "%s %d Fixed rank: %d\n", __FILE__, __LINE__, global_fixed_rank);
    }

    if(0 && num_mpi_ranks == 1 && initial_maxrank > N){ //FIXME Enable for distributed mem
        fprintf(stderr,"%s %d %d\t|N:%d is less than actual maxrank:%d\n", __FILE__, __LINE__, MORSE_My_Mpi_Rank(), N, initial_maxrank);
        exit(1);
    }
    int set_diag = 0;

    /* Save A for check */
    PROGRESS("pasting original dense descAD into Adense and Adense2 started");
    // Adense: original dense problem.
    PASTE_TILE_TO_LAPACK( descDense, Adense, check, double, LDA, M );
    double one = 1.0, zero = 0.0, minusone = -1.0, diagVal = M;
    double* swork = NULL;
    //double* cp_L_Adense =  calloc(LDA*M, sizeof(double)); 
    if(check){
        swork  = calloc(2*M, sizeof(double));
        {size_t i, j;
            double* orgAdense = calloc(LDA*M, sizeof(double));
            for(j = 0; j < M; j++){
                for(i = 0; i < M; i++){
                    orgAdense[j*LDA+i] = Adense[j*LDA+i];
                }
            }
            int info = LAPACKE_dpotrf_work(
                    LAPACK_COL_MAJOR,
                    'L',
                    M, orgAdense, LDA);
            if(0 && info != 0){ //FIXME
                fprintf(stderr, "%s\t|%d\t|Error in LAPACK potrf. info:%d, This errors means "
                        "that the matrix generated is not positive definite\n", __FILE__, __LINE__, info);
            }
            for(j = 0; j < M; j++){
                for(i = 0; i < j; i++){
                    orgAdense[j*LDA+i] = zero;
                }
            }
            /*for(j = 0; j < M; j++) {                */
                /*for(i = 0; i < M; i++){*/
                    /*cp_L_Adense[j*LDA+i] = orgAdense[j*LDA+i];*/
                /*}*/
            /*}*/
            if(main_print_mat ){printf("L of Adense\n");printmat(orgAdense,M,M,LDA,MB, MB);}
            double normOrgAdense = 0.0;
            /*HICMA_znormest(M, M, orgAdense, &normOrgAdense, swork);*/
            /*printf("norm_L_OrgAdense:%e\n",normOrgAdense);*/
            free(orgAdense);
        }
    }
    PASTE_TILE_TO_LAPACK( descDense, Adense2, check, double, LDA, M );
    PROGRESS("pasting original dense descAD into Adense and Adense2 finished");
    PROGRESS("potrf started");
    START_TIMING();
    HICMA_zpotrf_Tile(MorseLower, descAUV, descAD, descArk, fixedrank, maxrank, fixedacc );
    STOP_TIMING();
    fflush(stderr);
    fflush(stdout);
    PROGRESS("potrf finished");
    if(check){
        HICMA_zuncompress(MorseLower, descAUV, descDense, descArk);
        HICMA_zdiag_vec2mat(descAD, descDense);
        PASTE_CODE_FREE_MATRIX( descAD ); //@KADIRLBL001  
        descAD = descDense; // descAD was only diagonals.
        // After this line, descAD is dense matrix containing approximate L
        // So no need to adapt below code for descAD containg only diagonals.
    }
    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descArk, Ark_final, 1, double, MT, NT );
        if(MORSE_My_Mpi_Rank()==0){
            sprintf(rankfile, "%s-2", rankfile);
            fwrite_array(descArk->m, descArk->n, descArk->m, Ark_final, rankfile);
            HICMA_stat_t hicma_statrk_final;
            zget_stat(MorseLower, Ark_final, MT, NT, MT,  &hicma_statrk_final);
            printf("final_ranks:");
            zprint_stat(hicma_statrk_final);
            fflush(stderr);
            fflush(stdout);
        }
    }

    int check_dense = 0;
    int check_app = 1;
    if(check == 0){
        check_dense = check_app = 0;
    }
    if(check_app ) {
        PROGRESS("checking accuracy");
        if( MORSE_My_Mpi_Rank()==0){
#ifndef COMPLEX
            if(main_print_mat){printf("Adense2\n");printmat(Adense2,M,M,LDA,MB, MB);}
            double normA;
            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = 0; i < j; i++){
                        Adense2[j*LDA+i] = zero;
                    }
                }
            }
            PROGRESS("normaA started");
            HICMA_znormest(M, M, Adense2, &normA, swork);
            // Ahicma: result of TLR potrf
            PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, double, LDA, M );
            /*if(0){size_t i,j;*/
                /*for(j = 0; j < M; j++) {                */
                    /*for(i = 0; i < M; i++){*/
                        /*Ahicma[j*LDA+i] = cp_L_Adense[j*LDA+i];*/
                    /*}*/
                /*}*/
            /*}*/
            double normAhicma = 0.0;
            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = 0; i < j; i++){
                        Ahicma[j*LDA+i] = zero;
                    }
                }
                double* orgAhicma = calloc(LDA*M, sizeof(double));
                for(j = 0; j < M; j++){
                    for(i = 0; i < M; i++){
                        orgAhicma[j*LDA+i] = Ahicma[j*LDA+i];
                    }
                }
                HICMA_znormest(M, M, orgAhicma, &normAhicma, swork);
                free(orgAhicma);
            }
            if(set_diag){size_t j; for(j = 0; j < M; j++){ Ahicma[j*LDA+j] = diagVal; } }
            if(main_print_mat){printf("Ahicma\n");printmat(Ahicma,M,M,LDA, MB, MB);}
            //LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', M, Ahicma, LDA);
            // AhicmaT: transpose of Ahicma
            PROGRESS("copy descAd into AhicmaT started");
            PASTE_TILE_TO_LAPACK( descAD,  AhicmaT, check, double, LDA, M );

            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = 0; i < j; i++){
                        Adense[j*LDA+i] = zero;
                    }
                }
            }

            if(main_print_mat){printf("Ahicma-upperzero\n");printmat(Ahicma,M,M,LDA, MB, MB);}
            PROGRESS("Transpose A started");
            LAPACKE_dge_trans(LAPACK_COL_MAJOR, M, M, Ahicma, LDA, AhicmaT, LDA);
            if(main_print_mat){printf("AhicmaT\n");printmat(AhicmaT,M,M,LDA, MB, MB);}
            PROGRESS("TRMM started");
            cblas_dtrmm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, M, one, Ahicma, LDA, AhicmaT, LDA);
            if(main_print_mat){printf("Ahicma*AhicmaT\n");printmat(AhicmaT,M,M,LDA, MB, MB);}
            //double tmpnorm;normest(M, M, AhicmaT, &tmpnorm, swork);printf("tmpnorm:%e\n",tmpnorm);
            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = 0; i < j; i++){
                        AhicmaT[j*LDA+i] = zero;
                    }
                }
            }

            size_t nelm = M * M;
            if(main_print_mat)printf("nelm:%zu M:%d N:%d\n", nelm, M, N);
            PROGRESS("DAXPY started");
            cblas_daxpy(nelm, minusone, AhicmaT, 1, Adense, 1);
            if(main_print_mat){printf("Adense-(Ahicma*AhicmaT)\n");printmat(Adense,M,M,LDA, MB, MB);}

            double normDenseAppDiff;
            PROGRESS("Norm of difference started");
            HICMA_znormest(M, M, Adense, &normDenseAppDiff, swork);
            double accuracyDenseAppDiff = normDenseAppDiff/normA;
            //printf("normA:%.2e normDenseAppdiff:%.2e Accuracy: %.2e\n", normA, normDenseAppDiff,  accuracyDenseAppDiff);
            dparam[IPARAM_RES] = normDenseAppDiff;
            dparam[IPARAM_ANORM] = normA;
            dparam[IPARAM_XNORM] = normA;
            dparam[IPARAM_BNORM] = normAhicma;
#endif
        } else {
            PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, double, LDA, M );
            PASTE_TILE_TO_LAPACK( descAD,  AhicmaT, check, double, LDA, M );
        }
        PROGRESS("checking accuracy is finished");
    }

    PASTE_CODE_FREE_MATRIX( descAUV );
    PROGRESS("descAUV is freed");
    if(check == 0) { // If there is no check, then descAD and descDense are different. Refer to @KADIRLBL001
        PASTE_CODE_FREE_MATRIX( descAD );
        PROGRESS("descAD is freed");
    }
    PASTE_CODE_FREE_MATRIX( descArk );
    PROGRESS("descArk is freed");
    PASTE_CODE_FREE_MATRIX( descDense );
    PROGRESS("descDense is freed");
    PROGRESS("freed descs");
    return 0;
}
