/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file time_dgemm_tile.c
 *
 * This file shows how to generate tile low-rank (TLR) matrices and perform matrix multiplication.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2019-11-21
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

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include "timing.h"
#include <hicma_constants.h>
#include <hicma_struct.h>
#include <include/hicma_d.h>

#include "timing_auxiliary.h"

#include "starpu.h"

#include <assert.h>
#include "misc/auxdescutil.h"

#include <hicma.h>
#include <hicma_common.h>
#include "timing_dauxiliary.h"

#include "flop_util_structs.h"
extern flop_counter counters[FLOP_NUMTHREADS];
#include "flop_util.h"

//#if __has_include("starpu_mpi.h")
#ifdef HICMA_USE_MPI
#define __ENABLE_MPI
#include "starpu_mpi.h"
#endif

#include <math.h>
#include <time.h>

#undef  CBLAS_SADDR
#define CBLAS_SADDR(_val) (_val)

int print_progress = 1;   // Print progress about the execution
char datebuf[128];
time_t timer;
struct tm* tm_info;
#undef PROGRESS
#define PROGRESS(str)

int store_only_diagonal_tiles = 0;
int global_check = 0;
extern int global_always_fixed_rank;
extern int global_fixed_rank;
int num_mpi_ranks;
int diag_nrows = 0;
int main_print_index = 0;
extern int print_index;
int main_print_mat =0;
extern int print_mat;
int use_scratch = 1; // Use scratch memory provided by starpu
int calc_rank_stat = 1; 
extern int reorder_inner_products;
int comp_idrank (const void * elem1, const void * elem2) 
{
    idrank f = *((idrank*)elem1);
    idrank s = *((idrank*)elem2);
    if (f.rank > s.rank) return  1;
    if (f.rank < s.rank) return -1;
    return 0;
}
extern idrank*** iporder; //inner product order

double timediff(struct timeval begin, struct timeval end){
    double elapsed = (end.tv_sec - begin.tv_sec) +
        ((end.tv_usec - begin.tv_usec)/1000000.0);
    return elapsed;
}
// This function isused to compare to descriptors in terms of size and numerical values.
// Do not use this function with MPI
// FIXME: This function will be moved to aux/ folder
void check_same(HICMA_desc_t *descL, HICMA_desc_t *descR, char diag, char lower_triangle){
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
RunTest(int *iparam, double *dparam, hicma_time_t *t_, char* rankfile)
{
    // print progress info only on ROOT process
    if(HICMA_My_Mpi_Rank() != 0)
        print_progress = 0;
    PROGRESS("RunTest started");
    // set alpha and beta which will be used in gemm
    double alpha = 1.0, beta = 1.0;
    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 0;

    //HICMA_set_use_fast_hcore_dgemm();

    // Gemm requires more descriptors
    //chameleon/runtime/starpu/control/runtime_descriptor.c
    HICMA_user_tag_size(31,26);
    //HICMA_user_tag_size(31,28);

    global_always_fixed_rank = iparam[IPARAM_HICMA_ALWAYS_FIXED_RANK];

    // get parameters coming from command line
    PASTE_CODE_IPARAM_LOCALS( iparam );
    // set global variable so that p.. files can fill dense matrix
    global_check = check;
    // calculate total number of mpi processes
    num_mpi_ranks = P*Q;
    print_index = iparam[IPARAM_HICMA_PRINTINDEX];
    print_mat   = iparam[IPARAM_HICMA_PRINTMAT];
    int64_t _nb = iparam[IPARAM_NB];
    LDB = hicma_max(K, iparam[IPARAM_LDB]);
    LDC = hicma_max(M, iparam[IPARAM_LDC]);
    int hicma_maxrank  = iparam[IPARAM_HICMA_MAXRANK];
    reorder_inner_products = iparam[IPARAM_HICMA_REORDER_INNER_PRODUCTS];

    unsigned long*    opcounters  = NULL; // count all operations performed by a core
    int nelm_opcounters = num_mpi_ranks*thrdnbr;
    opcounters     = calloc(nelm_opcounters,     sizeof(unsigned long)); //TODO free
    flop_util_init_counters(thrdnbr);

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
    // TODO: since HICMA_dgytlr_Tile does not handle NULL dense pointer, have to create them, better condition it with 'check' variable
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descA, check, double, HicmaRealDouble, LDA, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descB, check, double, HicmaRealDouble, LDB, M, M );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descC, check, double, HicmaRealDouble, LDC, M, M );
    PROGRESS("desc..'s are allocated");

    // gemm will not have dense diagonals in RND TILED. So remove desc..D's in the future. 
    // since HICMA_dgytlr_Tile does not handle NULL dense pointer, have to create them, better condition it with 'check' variable
    /*PASTE_CODE_ALLOCATE_MATRIX_TILE( descAD, 1, double, HicmaRealDouble, LDA, M, M );*/
    /*PASTE_CODE_ALLOCATE_MATRIX_TILE( descBD, 1, double, HicmaRealDouble, LDB, M, M );*/
    /*PASTE_CODE_ALLOCATE_MATRIX_TILE( descCD, 1, double, HicmaRealDouble, LDC, M, M );*/
    NB = saveNB;


    //U and V's
    saveN = N;
    N = N * 2;
    saveNB = NB;
    NB = NB * 2;
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUV, 1, double, HicmaRealDouble, LDC, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descBUV, 1, double, HicmaRealDouble, LDC, M, N );
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descCUV, 1, double, HicmaRealDouble, LDC, M, N );
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
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descArk, 1, double, HicmaRealDouble, ldrk, mt, mt);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descBrk, 1, double, HicmaRealDouble, ldrk, mt, mt);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descCrk, 1, double, HicmaRealDouble, ldrk, mt, mt);
    MB = bigMB;
    NB = bigNB;
    PROGRESS("desc..rk's are allocated");

    int _k = iparam[IPARAM_RK]; //genargs->k
    //double _acc = pow(10, -1.0*iparam[IPARAM_ACC]);
    double _acc = dparam[IPARAM_HICMA_ACCURACY_THRESHOLD];

    char sym;
    sym = 'N';
    sym = 'S'; // FIXME
    int probtype = iparam[IPARAM_HICMA_STARSH_PROB];
    int maxrank  = iparam[IPARAM_HICMA_STARSH_MAXRANK];
    double ddecay = dparam[IPARAM_HICMA_STARSH_DECAY];
    int  initial_maxrank, final_maxrank;
    double initial_avgrank, final_avgrank;
    HICMA_problem_t hicma_problem;

    hicma_problem.ndim = 2;
    //BEGIN: rndtiled
    if(iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_RND) {
        hicma_problem.noise    = 0.0;
    }
    //END:   rndtiled

    //BEGIN: geostat
    double theta[3] = {1, 0.1, 0.5};
    if(iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_GEOSTAT) {
        hicma_problem.theta = theta;
        hicma_problem.noise = 0.0;
    }
    //END: geostat

    //BEGIN: ss
    if(iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_SS) {
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
    //BEGIN: st-3D-exp
    if(iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_ST_3D_EXP) {
        /*
        // from lorapo
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        double sigma = 1.0;
        int ndim = 3;
        kernel = starsh_ssdata_block_exp_kernel_3d; 
        info = starsh_ssdata_generate((STARSH_ssdata **)&data, N, ndim,
                beta, nu, noise,
                place, sigma);
         */
        // Correlation length
        hicma_problem.beta  = 0.1;
        //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

        // Smoothing parameter for Matern kernel
        hicma_problem.nu    = 0.5;
        // Shift added to diagonal elements
        hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
        hicma_problem.noise = 5.e-4; //works for 640K
    }
    //END: st-3D-exp
    //BEGIN: edsin
    if(iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_EDSIN) {
        // Wave number, >= 0
        hicma_problem.wave_k = dparam[IPARAM_HICMA_STARSH_WAVE_K];
        hicma_problem.diag = 0; 
    //printf("%s %d: %g\n", __FILE__, __LINE__, hicma_problem.wave_k);
    }
    //END: edsin
    PROGRESS("generating coordinates started");
    generate_problem(probtype, sym, ddecay, M, MB, mt, nt, &hicma_problem);
    PROGRESS("generating coordinates ended");

	// enforce compression of diagonal tiles
	int compress_diag = 1;
    PROGRESS("dgytlr starting");
    //desc[A,B,C] are original problem if global check is enabled
    //desc[A,B,C]D are original problem. Temporary buffers should be used in hcore_dgytlr TODO Do not use desc..D in this gemm code
    /*HICMA_dgytlr_Tile(HicmaUpperLower, descAUV, descAD, descArk, 0, maxrank, _acc, compress_diag,  descA);*/
    /*HICMA_dgytlr_Tile(HicmaUpperLower, descBUV, descBD, descBrk, 0, maxrank, _acc, compress_diag,  descB);*/
    /*HICMA_dgytlr_Tile(HicmaUpperLower, descCUV, descCD, descCrk, 0, maxrank, _acc, compress_diag,  descC);*/
   
    if(check == 1){
        HICMA_dhagdm_Tile(HicmaUpperLower, descA);
        HICMA_dhagdm_Tile(HicmaUpperLower, descB);
        HICMA_dhagdm_Tile(HicmaUpperLower, descC);
    }

    HICMA_dhagcm_Tile(HicmaUpperLower, descAUV, descArk, M, M, MB,  MB, maxrank, _acc); //FIXME sizes
    HICMA_dhagcm_Tile(HicmaUpperLower, descBUV, descBrk, M, M, MB,  MB, maxrank, _acc); //FIXME sizes
    HICMA_dhagcm_Tile(HicmaUpperLower, descCUV, descCrk, M, M, MB,  MB, maxrank, _acc); //FIXME sizes

    PROGRESS("dgytlr finished");


    double* _Ark = NULL;
    double* _Brk = NULL;
    if(calc_rank_stat == 1 || reorder_inner_products == 1){
        PASTE_TILE_TO_LAPACK( descArk, Ark, 1, double, MT, NT );
        PASTE_TILE_TO_LAPACK( descBrk, Brk, 1, double, MT, NT );
        _Ark = Ark;
        _Brk = Brk;
    }
    if(reorder_inner_products == 1){
        iporder = calloc(MT, sizeof(idrank**));
        for(int i=0; i < MT; i++){
            iporder[i] = calloc(NT, sizeof(idrank*));
            for(int j=0; j < NT; j++){
                iporder[i][j] = calloc(MT, sizeof(idrank));
                for(int k=0; k < MT; k++){ // FIXME if matrices are rectangular
                    int rankA = _Ark[k*ldrk+i];
                    int rankB = _Brk[j*ldrk+k];
                    iporder[i][j][k].id = k;
                    iporder[i][j][k].rank = rankA > rankB ? rankA : rankB;
                    //printf("%d,%d:%d\t%d,%d:%d [%d,%d]\n", i, k, rankA, k, j, rankB, iporder[i][j][k].id, iporder[i][j][k].rank);
                }
				/** sort tiles in ascending order of ranks */
                qsort(iporder[i][j], MT, sizeof(idrank), comp_idrank);
                for(int _k=0; _k < MT; _k++){ // FIXME if matrices are rectangular
					int k = iporder[i][j][_k].id;
                    int rankA = _Ark[k*ldrk+i];
                    int rankB = _Brk[j*ldrk+k];
                    //printf("%d,%d:%d\t%d,%d:%d\n", i, k, rankA, k, j, rankB);
                }
				//getc(stdin);
            }
        }
    }
    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descArk, Ark_initial, 1, double, MT, NT );
        if(HICMA_My_Mpi_Rank()==0){

            /*sprintf(rankfile, "%s-1", rankfile);*/
            //fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);
            print_array(descArk->m, descArk->n, descArk->m, Ark_initial, stdout);

            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaUpperLower, Ark_initial, mt, mt, mt,  &hicma_statrk_initial);
            printf("Ainitial_ranks:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }
    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descCrk, Crk_initial, 1, double, MT, NT );
        if(HICMA_My_Mpi_Rank()==0){

            /*sprintf(rankfile, "%s-1", rankfile);*/
            /*fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);*/

            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaUpperLower, Crk_initial, mt, mt, mt,  &hicma_statrk_initial);
            printf("Cinitial_ranks:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }

    int set_diag = 0;
    double diagVal = M;

    /* Save A for check */
    PROGRESS("pasting original dense descC into C2 started");
    // Adense: original dense problem.
    PASTE_TILE_TO_LAPACK( descC, C2, check, double, LDC, M );
    PROGRESS("pasting original dense descC into C2 finished");

    PROGRESS("gemm started");
    START_TIMING();
    HICMA_dgemm_Tile( HicmaNoTrans, HicmaNoTrans, alpha,  //TODO
            descAUV, descArk,
            descBUV, descBrk,
            beta,
            descCUV, descCrk, _k, maxrank, _acc
            );
    STOP_TIMING();
    fflush(stderr);
    fflush(stdout);
    PROGRESS("gemm finished");

    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descCrk, Crk_final, 1, double, MT, NT );
        if(HICMA_My_Mpi_Rank()==0){

            /*sprintf(rankfile, "%s-1", rankfile);*/
            /*fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);*/

            HICMA_stat_t hicma_statrk_final;
            dget_stat(HicmaUpperLower, Crk_final, mt, mt, mt,  &hicma_statrk_final);
            printf("Cfinal_ranks:");
            dprint_stat(hicma_statrk_final);
            fflush(stderr);
            fflush(stdout);
        }
    }
    assert(thrdnbr < FLOP_NUMTHREADS);
    int myrank = HICMA_My_Mpi_Rank();
    for(int i = 0; i < thrdnbr; i++){
        flop_counter res = counters[i];
        unsigned long totflop = res.potrf+res.trsm+res.syrk+res.update;
        //printf("myrank:%d thread:%d %lu\n", myrank, i, totflop);
        opcounters[myrank * thrdnbr + i] = totflop;
    }

    unsigned long* allopcounters = opcounters;
#ifdef __ENABLE_MPI
    allopcounters     = calloc(nelm_opcounters,     sizeof(unsigned long)); //TODO free
    MPI_Reduce(opcounters, allopcounters, nelm_opcounters, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);  
#endif
    unsigned long totflop = 0;
    if(HICMA_My_Mpi_Rank() == 0) {
        unsigned long sum_thread = 0;
        printf("nop_thread %d %d\n", num_mpi_ranks, thrdnbr);
        for(int i = 0; i < num_mpi_ranks; i++){
            for(int j = 0; j < thrdnbr; j++){
                printf("%lu ", allopcounters[i*thrdnbr+j]);
                sum_thread += allopcounters[i*thrdnbr+j];
            }
            printf("\n");
        }
        totflop += sum_thread;
    }
    if(HICMA_My_Mpi_Rank() == 0)
    {
        /** prints number of flops */
        char str[1024];
        str[1023] = '\0';
        //printf("\t\tPOTRF\tTRSM \tSYRK\tGEMM\t\tTOTFLOP\t\tTOTGFLOP\tGFLOP/S\t\tTIME(s)\n");
        printf("\t\tTOTFLOP\t\tTOTGFLOP\tGFLOP/S\t\tTIME(s)\n");
        printf("ReShg\t");
        printf("%lu\t", totflop);
        double totgflop = totflop/(1024.0*1024*1024);
        printf("%g\t", totgflop);
        double totgflops = totgflop/t;
        printf("%g\t", totgflops);
        printf("%g", t);
        printf("\n");
    }

    if(check){

        PASTE_TILE_TO_LAPACK( descA, A, check, double, LDA, M );
        PASTE_TILE_TO_LAPACK( descB, B, check, double, LDB, M );
        HICMA_duncompress(HicmaUpperLower, descCUV, descC, descCrk);

        PASTE_TILE_TO_LAPACK( descC, C, check, double, LDC, M );


        dparam[IPARAM_RES] = hicma_d_check_gemm( HicmaNoTrans, HicmaNoTrans, M, M, M,
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

    PASTE_CODE_FREE_MATRIX( descArk );
    PASTE_CODE_FREE_MATRIX( descBrk );
    PASTE_CODE_FREE_MATRIX( descCrk );
    PROGRESS("desc..rk's are freed");

    if(check == 1) {
        PASTE_CODE_FREE_MATRIX( descA );
        PASTE_CODE_FREE_MATRIX( descB );
        PASTE_CODE_FREE_MATRIX( descC );
    }
    PROGRESS("desc..'s are freed");
    PROGRESS("freed descs");
    return 0;
}
