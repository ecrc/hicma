/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file time_dpotrf_tile.c
 *
 * This file shows how to generate tile low-rank (TLR) matrix and factorize it using Cholesky factorization.
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
 * The meaning of the descriptors:
 * - AUV: U and V, side by side
 * - AD : U*V
 * - A  : the original, non-approximated problem
 * - Ark: rank of U and V, each tile of the matrix is a single integer in fact a double.
 *
 **/

#include <stdio.h>
#include <time.h>
#include <sys/time.h>

#include <hicma.h>
#include "timing.h"
#include <hicma_constants.h>
#include <hicma_struct.h>
#include <include/hicma_d.h>
#include "timing_auxiliary.h"
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
#include "starsh-rbf.h"

#include <assert.h>
#include "hicma_common.h"
#include "misc/auxcompute_d.h"
#include "misc/auxdescutil.h"
#include "hicma.h"

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

extern double rad; // RBF scaling factor
extern double reg; // RBF regularization value
extern double denst; //RBF density
extern char *mesh_file;
extern char *interpl_file;

int print_progress = 0;   // Print progress about the execution
char datebuf[128];
time_t timer;
struct tm *tm_info;
#define PROGRESS(str) \
    if(print_progress){ \
        int myrank = HICMA_My_Mpi_Rank();\
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
extern int global_always_fixed_rank;
extern int global_fixed_rank;
int global_omit_computation = 1;
int num_mpi_ranks;
int run_potrf = 1;
int diag_nrows = 0;
int main_print_index = 0;
extern int print_index;
extern int print_index_end;
int main_print_mat = 0;
extern int print_mat;
int use_scratch = 1; // Use scratch memory provided by starpu
int calc_rank_stat = 0;
double total_compute_time = 0.0, tcompute = 0.0;
double total_compress_time = 0.0, tcompress = 0.0;

double timediff(struct timeval begin, struct timeval end) {
    double elapsed = (end.tv_sec - begin.tv_sec) +
                     ((end.tv_usec - begin.tv_usec) / 1000000.0);
    return elapsed;
}

/*This is for fine timing for computation only or compression 
  only over several iterations */
#define START_TIMING_COMPUTE()                \
  tcompute = -HICMA_RUNTIME_get_time();

#define STOP_TIMING_COMPUTE()                 \
  tcompute += HICMA_RUNTIME_get_time();            \
  total_compute_time = total_compute_time + tcompute;


#define START_TIMING_COMPRESSION()                \
  tcompress = -HICMA_RUNTIME_get_time();

#define STOP_TIMING_COMPRESSION()                 \
  tcompress += HICMA_RUNTIME_get_time();            \
  total_compress_time = total_compress_time + tcompress;


int
RunTest(int *iparam, double *dparam, hicma_time_t *t_, char *rankfile) {
    // print progress info only on ROOT process
    if (HICMA_My_Mpi_Rank() != 0)
        print_progress = 0;
    PROGRESS("RunTest started");

    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 1;
    //chameleon/runtime/starpu/control/runtime_descriptor.c
    //HICMA_user_tag_size(31,26);
    //HICMA_user_tag_size(31,29);
    // HICMA_user_tag_size(31,27);// When I added tile_to_lapack for descArk, I got not enough number of desc error
    HICMA_user_tag_size(64, 50);
    // get parameters coming from command line
    PASTE_CODE_IPARAM_LOCALS(iparam);
    // set global variable so that p.. files can fill dense matrix
    global_check = check;
    // calculate total number of mpi processes (it is not used for now)
    num_mpi_ranks = P * Q;

    //This is for batch mode because each mpi node will loop over small number of batches and compute TLR cholesky, so all allocation will be the same across all nodes
    P = 1;
    Q = 1;

    print_index = iparam[IPARAM_HICMA_PRINTINDEX];
    print_index_end = iparam[IPARAM_HICMA_PRINTINDEXEND];
    print_mat = iparam[IPARAM_HICMA_PRINTMAT];
    int64_t _nb = iparam[IPARAM_NB];
    LDA = hicma_max(M, iparam[IPARAM_LDA]);
    int hicma_maxrank = iparam[IPARAM_HICMA_MAXRANK];
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
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descAD, 1, double, HicmaRealDouble, LDA, M, N);
    N = saveN;
    P = saveP;
    Q = saveQ;
    PROGRESS("descAD is allocated");

    size_t ncols_Dense;
    size_t ld_Dense;
    int saveMB = MB;
    if (check == 0) {
        ncols_Dense = MT;
        MB = NB = 1;
        ld_Dense = MT;
    } else {
        ncols_Dense = M;
        ld_Dense = M;
    }
    /*descDense is full matrix if numerical accuracy will be checked.
     * Otherwise it is MB-by-MB matrix with 1-by-1 tiles */
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descDense, 1, double, HicmaRealDouble, ld_Dense, ncols_Dense, ncols_Dense);

    MB = saveMB;
    NB = 3;
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descB, solve, double, HicmaRealDouble, M, M, 3);
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descBcpy, check_solve, double, HicmaRealDouble, M, M, 3);

    MB = saveMB;
    NB = MB;
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descmat, check_solve, double, HicmaRealDouble, M, M, M);
    if (check == 0) {
        MB = saveMB;
    } else {
    }
    NB = saveNB;

    PROGRESS("descDense is allocated");

    int MTMB = MT * MB; // roundup number of rows/columns for AUV 
    int nrows_AUV = MTMB;
    int ld_AUV = MTMB;
    // allocate descUV
    saveN = N;
    N = N * 2;
    saveNB = NB;
    NB = NB * 2;
    //printf("N:%d NB:%d\n", N, NB);
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descAUV, 1, double, HicmaRealDouble, ld_AUV, nrows_AUV, N);

    N = saveN;
    NB = saveNB;
    PROGRESS("descAUV is allocated");

    /* tile dimension of rank descriptor must be 1 */
    /* when LD for rk matrices is 1, program exits*/
    int bigMB = MB;
    int bigNB = NB;
    MB = NB = 1;
    PASTE_CODE_ALLOCATE_MATRIX_TILE(descArk, 1, double, HicmaRealDouble, MT, MT, NT);
    PROGRESS("descA's are allocated");
    MB = bigMB;
    NB = bigNB;

    int diag_dense = 1;
    int fixedrank = iparam[IPARAM_RK]; //genargs->k
    //double fixedacc = pow(10, -1.0*iparam[IPARAM_ACC]);
    double fixedacc = dparam[IPARAM_HICMA_ACCURACY_THRESHOLD];

    char sym;
    if (run_potrf)
        sym = 'S';
    else
        sym = 'N';
    int probtype = iparam[IPARAM_HICMA_STARSH_PROB];
    int maxrank = iparam[IPARAM_HICMA_STARSH_MAXRANK];
    //double ddecay = pow(10, -1.0*iparam[IPARAM_HICMA_STARSH_DECAY]);
    double ddecay = dparam[IPARAM_HICMA_STARSH_DECAY];
    START_TIMING();
    /*
       Each MPI node will start from it is mpi rank for the first file to read. Then,  
       each MPI node will follow 1D cyclic distrbution when processing next mesh file
     */

    //compute total number of batches
    double batches = iparam[IPARAM_NUMOBJ] / iparam[IPARAM_NUMSUBOBJ];
    //compute number of batches per node
    int batches_per_node = (int) ceil(batches / (num_mpi_ranks));


    //TODO indentation, end of curly barckets
    int filenumber = HICMA_My_Mpi_Rank();
    for (int k = 0; k < batches_per_node; k++) {

        int initial_maxrank, final_maxrank;
        double initial_avgrank, final_avgrank;
        HICMA_problem_t hicma_problem;

        hicma_problem.ndim = 2;

        //BEGIN: rndtiled
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_RND) {
            hicma_problem.noise = 1.0; //value added to diagonal
        }
        //END:   rndtiled

        //BEGIN: geostat
        //double theta[3] = {1, 0.1, 0.5}; //initially
        double theta[3] = {
                1.0, //sigma
                0.01, //beta
                10.0 //nu Aleks used 10.0 in his paper
        };
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_GEOSTAT) {
            hicma_problem.theta = theta;
            hicma_problem.noise = 0.0;
            hicma_problem.noise = 1.e-2;
            hicma_problem.kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
        }
        //END: geostat

        //BEGIN: ss
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_SS) {
            //sigma=1.0 default value line 193 of stars-h/src/applications/spatial.c 
            // Correlation length
            hicma_problem.beta = 0.1;
            //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

            // Smoothing parameter for Matern kernel
            hicma_problem.nu = 0.5;

            // Shift added to diagonal elements
            hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
            hicma_problem.noise = 5.e-4; //works for 640K but does not work for 10M
            hicma_problem.noise = 1.e-2; //
        }
        //END: ss
        //BEGIN: st-3D-exp
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_ST_3D_EXP) {
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
            hicma_problem.beta = 0.1;
            //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

            // Smoothing parameter for Matern kernel
            hicma_problem.nu = 0.5;
            // Shift added to diagonal elements
            hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
            hicma_problem.noise = 5.e-4; //works for 640K
            hicma_problem.noise = 1.e-2; //
        }
        //END: st-3D-exp
        //BEGIN: st-3D-sqexp
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_ST_3D_SQEXP) {
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
            hicma_problem.beta = 0.1;
            //If fixed rank is required set beta=1 and a sample case will be like this nb=25 maxrank=10 m=2500 So ranks will decrease.

            // Smoothing parameter for Matern kernel
            hicma_problem.nu = 0.5;
            // Shift added to diagonal elements
            hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
            hicma_problem.noise = 5.e-4; //works for 640K
            hicma_problem.noise = 1.e-2; //
        }
        //END: st-3D-exp
        //BEGIN: edsin
        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_EDSIN) {
            // Wave number, >= 0
            hicma_problem.wave_k = dparam[IPARAM_HICMA_STARSH_WAVE_K];
            hicma_problem.diag = M;
            //printf("%s %d: %g\n", __FILE__, __LINE__, hicma_problem.wave_k);
        } //END: edsin

        //RBF Unstructured Mesh Deformation for 3D problems

        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_3D_RBF_VIRUS) {
            hicma_problem.kernel_type = iparam[IPARAM_RBFKERNEL]; // RBF kernel_type
            hicma_problem.reg = 1 + fixedacc * 10;  //RBF regularization value
            hicma_problem.isreg = 1;
            hicma_problem.rad = rad;  // RBF scaling factor
            hicma_problem.denst = denst;
            hicma_problem.mesh_points = M;
            hicma_problem.mordering = iparam[IPARAM_ORDER];
            hicma_problem.numobj = iparam[IPARAM_NUMSUBOBJ];  // how many subobjects (e.g. number of subviruses in the batch)

            /* In batch mode user has to define only path to directory that containes mesh files
               File name has to be in sequence starting from 0.txt, 1.txt, and so on
               The following line will construct the whole path to file including file name
             */
            char filenum[100];
            char filepath[1000];
            sprintf(filenum, "%d", filenumber);
            strcpy(filepath, mesh_file);
            strcat(filepath, filenum);
            strcat(filepath, ".txt");
            hicma_problem.mesh_file = filepath;   // path to mesh file
            if (0)
                printf("\nI am mpi:%d start batch:%d and batches per node:%d and total batches:%f\n",
                       HICMA_My_Mpi_Rank(), k, batches_per_node, batches);
            if (0)printf("\nI am mpi:%d, my file is:%s\n", HICMA_My_Mpi_Rank(), filepath);
        }


        PROGRESS("generating coordinates started");
        struct timeval tvalBefore, tvalAfter;  // removed comma
        gettimeofday(&tvalBefore, NULL);
        generate_problem(probtype, sym, ddecay, M, MB, MT, NT, &hicma_problem);
        gettimeofday(&tvalAfter, NULL);
        if (HICMA_My_Mpi_Rank() == 0) {
            printf("Tproblem:%g\n",
                   (tvalAfter.tv_sec - tvalBefore.tv_sec)
                   + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0
            );
            fflush(stderr);
            fflush(stdout);
        }
        PROGRESS("generating coordinates ended");

        // DO NOT enforce compression of diagonal tiles
        int compress_diag = 0;

        START_TIMING_COMPRESSION();
        PROGRESS("nompi dgytlr starting");
        //descDense original problem
        gettimeofday(&tvalBefore, NULL);
        HICMA_dgytlr_Tile(HicmaLower, descAUV, descAD, descArk, 0, maxrank, fixedacc, compress_diag, descDense);
        gettimeofday(&tvalAfter, NULL);
        if (HICMA_My_Mpi_Rank() == 0) {
            printf("Tcompress:%g\n",
                   (tvalAfter.tv_sec - tvalBefore.tv_sec)
                   + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0
            );
            fflush(stderr);
            fflush(stdout);
        }
        PROGRESS("nompi dgytlr finished");
        fflush(stderr);
        fflush(stdout);
        /*return 0; //TODO*/
        STOP_TIMING_COMPRESSION();

        if (calc_rank_stat == 1) {
            PASTE_TILE_TO_LAPACK(descArk, Ark_initial, 1, double, MT, NT);
            if (HICMA_My_Mpi_Rank() == 0) {

                if (0) sprintf(rankfile, "%s_initialranks", rankfile);
                if (0) fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);
                print_array(descArk->m, descArk->n, descArk->m, Ark_initial, stdout);

                HICMA_stat_t hicma_statrk_initial;
                dget_stat(HicmaLower, Ark_initial, MT, NT, MT, &hicma_statrk_initial);
                printf("initial_ranks:");
                dprint_stat(hicma_statrk_initial);
                fflush(stderr);
                fflush(stdout);
            }
        }

        if (global_always_fixed_rank == 1) {
            fprintf(stderr, "%s %d Fixed rank: %d\n", __FILE__, __LINE__, global_fixed_rank);
        }

        if (0 && num_mpi_ranks == 1 && initial_maxrank > N) { //FIXME Enable for distributed mem
            fprintf(stderr, "%s %d %d\t|N:%d is less than actual maxrank:%d\n", __FILE__, __LINE__, HICMA_My_Mpi_Rank(),
                    N, initial_maxrank);
            exit(1);
        }
        int set_diag = 0;

        /* Save A for check */
        PROGRESS("pasting original dense descAD into Adense and Adense2 started");
        // Adense: original dense problem.
        PASTE_TILE_TO_LAPACK(descDense, Adense, check, double, LDA, M);
        double one = 1.0, zero = 0.0, minusone = -1.0, diagVal = M;
        double *swork = NULL;
        //double* cp_L_Adense =  calloc(LDA*M, sizeof(double)); 
        if (check) {
            swork = calloc(2 * M, sizeof(double));
            {
                size_t i, j;
                double *orgAdense = calloc(LDA * M, sizeof(double));
                for (j = 0; j < M; j++) {
                    for (i = 0; i < M; i++) {
                        orgAdense[j * LDA + i] = Adense[j * LDA + i];
                    }
                }
                int info = LAPACKE_dpotrf_work(
                        LAPACK_COL_MAJOR,
                        'L',
                        M, orgAdense, LDA);
                if (0 && info != 0) { //FIXME
                    fprintf(stderr, "%s\t|%d\t|Error in LAPACK potrf. info:%d, This errors means "
                                    "that the matrix generated is not positive definite\n", __FILE__, __LINE__, info);
                }
                for (j = 0; j < M; j++) {
                    for (i = 0; i < j; i++) {
                        orgAdense[j * LDA + i] = zero;
                    }
                }
                /*for(j = 0; j < M; j++) {                */
                /*for(i = 0; i < M; i++){*/
                /*cp_L_Adense[j*LDA+i] = orgAdense[j*LDA+i];*/
                /*}*/
                /*}*/
                if (main_print_mat) {
                    printf("L of Adense\n");
                    printmat(orgAdense, M, M, LDA, MB, MB);
                }
                double normOrgAdense = 0.0;
                /*HICMA_znormest(M, M, orgAdense, &normOrgAdense, swork);*/
                /*printf("norm_L_OrgAdense:%e\n",normOrgAdense);*/
                free(orgAdense);
            }
        }
        PASTE_TILE_TO_LAPACK(descDense, Adense2, check, double, LDA, M);
        PROGRESS("pasting original dense descAD into Adense and Adense2 finished");

        unsigned long *opcounters = NULL; // count all operations performed by a core
        int nelm_opcounters = num_mpi_ranks * thrdnbr;
        opcounters = calloc(nelm_opcounters, sizeof(unsigned long)); //TODO free
        flop_util_init_counters(thrdnbr);

        PROGRESS("potrf started");
        START_TIMING_COMPUTE();
        HICMA_dpotrf_Tile(HicmaLower, descAUV, descAD, descArk, fixedrank, maxrank, fixedacc);
        STOP_TIMING_COMPUTE();
        fflush(stderr);
        fflush(stdout);
        PROGRESS("potrf finished");

        if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_3D_RBF_VIRUS) {
            if (solve) {
                // problem is 3D
                HICMA_dgenrhs_Tile(descB);

                PASTE_TILE_TO_LAPACK(descB, descB_array, 1, double, M, 3);
                double *descB_array_cpy = (double *) malloc(M * 3 * sizeof(double));

                if (check_solve) {
//                    HICMA_dlacpy_Tile(HicmaUpperLower, descB, descBcpy);
                    LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'X', M, 3, descB_array, M, descB_array_cpy, M);
                }

                HICMA_dtrsmd_Tile(HicmaLeft, HicmaLower, HicmaNoTrans, HicmaNonUnit, 1, descAUV, descAD, descArk, descB,
                                  maxrank);
                HICMA_dtrsmd_Tile(HicmaLeft, HicmaLower, HicmaTrans, HicmaNonUnit, 1, descAUV, descAD, descArk, descB,
                                  maxrank);

                if (check_solve) {

                    HICMA_dgenmat_Tile(descmat);
//                    dparam[IPARAM_IANORM] = HICMA_dlange_Tile(HicmaInfNorm, descmat);
//                    dparam[IPARAM_IBNORM] = HICMA_dlange_Tile(HicmaInfNorm, descBcpy);
//                    dparam[IPARAM_IXNORM] = HICMA_dlange_Tile(HicmaInfNorm, descB);

                    PASTE_TILE_TO_LAPACK(descmat, descmat_array, 1, double, M, M);
                    dparam[IPARAM_IANORM] = LAPACKE_dlange(LAPACK_COL_MAJOR, 'i', M, M, descmat_array, M);
                    dparam[IPARAM_IBNORM] = LAPACKE_dlange(LAPACK_COL_MAJOR, 'i', M, 3, descB_array_cpy, M);
                    dparam[IPARAM_IXNORM] = LAPACKE_dlange(LAPACK_COL_MAJOR, 'i', M, 3, descB_array, M);

                    double alpha = 1.0, beta = -1.0;

//                    HICMA_dgemm_Tile(HicmaNoTrans, HicmaNoTrans, alpha, descmat, descB, beta, descBcpy);
                    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M, 3, M, alpha,
                                descmat_array, M, descB_array, M, beta, descB_array_cpy, M);

//                    dparam[IPARAM_IRNORM] = HICMA_dlange_Tile(HicmaInfNorm, descBcpy);
                    dparam[IPARAM_IRNORM] = LAPACKE_dlange(LAPACK_COL_MAJOR, 'i', M, 3, descB_array_cpy, M);

                    dparam[IPARAM_IRES] = dparam[IPARAM_IRNORM] /
                                          (dparam[IPARAM_IANORM] * dparam[IPARAM_IXNORM] + dparam[IPARAM_IBNORM]);
                    free(descB_array_cpy);
                    free(descmat_array);
                }
                free(descB_array);

            }
        }
        assert(thrdnbr < FLOP_NUMTHREADS);
        int myrank = HICMA_My_Mpi_Rank();
        for (int i = 0; i < thrdnbr; i++) {
            flop_counter res = counters[i];
            unsigned long totflop = res.potrf + res.trsm + res.syrk + res.update;
            //printf("myrank:%d thread:%d %lu\n", myrank, i, totflop);
            if (0)
                printf("myrank:%d thread:%d po:%lu tr:%lu sy:%lu gm:%lu\n", myrank, i, res.potrf, res.trsm, res.syrk,
                       res.update);
            opcounters[myrank * thrdnbr + i] = totflop;
        }

        unsigned long *allopcounters = opcounters;
#ifdef __ENABLE_MPI
        allopcounters     = calloc(nelm_opcounters,     sizeof(unsigned long)); //TODO free
        MPI_Reduce(opcounters, allopcounters, nelm_opcounters, MPI_UNSIGNED_LONG, MPI_MAX, 0, MPI_COMM_WORLD);
#endif
        unsigned long totflop = 0;
        if (HICMA_My_Mpi_Rank() == 0) {
            unsigned long sum_thread = 0;
            if (0) printf("nop_thread %d %d\n", num_mpi_ranks, thrdnbr);
            for (int i = 0; i < num_mpi_ranks; i++) {
                for (int j = 0; j < thrdnbr; j++) {
                    if (0) printf("%lu ", allopcounters[i * thrdnbr + j]);
                    sum_thread += allopcounters[i * thrdnbr + j];
                }
                if (0)printf("\n");
            }
            totflop += sum_thread;
        }
        if (HICMA_My_Mpi_Rank() == 0) {
            /** prints number of flops */
            char str[1024];
            str[1023] = '\0';
            //printf("\t\tPOTRF\tTRSM \tSYRK\tGEMM\t\tTOTFLOP\t\tTOTGFLOP\tGFLOP/S\t\tTIME(s)\n");
            printf("\t\tTOTFLOP\t\tTOTGFLOP\tGFLOP/S\t\tTIME(s)\n");
            printf("ReShg\t");
            printf("%lu\t", totflop);
            double totgflop = totflop / (1024.0 * 1024 * 1024);
            printf("%g\t", totgflop);
            double totgflops = totgflop / tcompute;
            printf("%g\t", totgflops);
            printf("%g", tcompute);
            printf("\n");
        }

        if (check) {
            HICMA_duncompress(HicmaLower, descAUV, descDense, descArk);
            HICMA_ddiag_vec2mat(descAD, descDense);
            PASTE_CODE_FREE_MATRIX(descAD); //@KADIRLBL001
            descAD = descDense; // descAD was only diagonals.
            // After this line, descAD is dense matrix containing approximate L
            // So no need to adapt below code for descAD containg only diagonals.
        }
        if (calc_rank_stat == 1) {
            PASTE_TILE_TO_LAPACK(descArk, Ark_final, 1, double, MT, NT);
            if (HICMA_My_Mpi_Rank() == 0) {
                if (0) sprintf(rankfile, "%s_finalranks", rankfile);
                if (0) fwrite_array(descArk->m, descArk->n, descArk->m, Ark_final, rankfile);
                print_array(descArk->m, descArk->n, descArk->m, Ark_final, stdout);
                HICMA_stat_t hicma_statrk_final;
                dget_stat(HicmaLower, Ark_final, MT, NT, MT, &hicma_statrk_final);
                printf("final_ranks:");
                dprint_stat(hicma_statrk_final);
                fflush(stderr);
                fflush(stdout);
            }
        }

        int check_dense = 0;
        int check_app = 1;
        if (check == 0) {
            check_dense = check_app = 0;
        }
        if (check_app) {
            PROGRESS("checking accuracy");
            if (HICMA_My_Mpi_Rank() == 0) {
#ifndef COMPLEX
                if (main_print_mat) {
                    printf("Adense2\n");
                    printmat(Adense2, M, M, LDA, MB, MB);
                }
                double normA;
                {
                    size_t i, j;
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            Adense2[j * LDA + i] = zero;
                        }
                    }
                }
                PROGRESS("normaA started");
                HCORE_dnormest(M, M, Adense2, &normA, swork);
                // Ahicma: result of TLR potrf
                PASTE_TILE_TO_LAPACK(descAD, Ahicma, check, double, LDA, M);
                /*if(0){size_t i,j;*/
                /*for(j = 0; j < M; j++) {                */
                /*for(i = 0; i < M; i++){*/
                /*Ahicma[j*LDA+i] = cp_L_Adense[j*LDA+i];*/
                /*}*/
                /*}*/
                /*}*/
                double normAhicma = 0.0;
                {
                    size_t i, j;
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            Ahicma[j * LDA + i] = zero;
                        }
                    }
                    double *orgAhicma = calloc(LDA * M, sizeof(double));
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < M; i++) {
                            orgAhicma[j * LDA + i] = Ahicma[j * LDA + i];
                        }
                    }
                    HCORE_dnormest(M, M, orgAhicma, &normAhicma, swork);
                    free(orgAhicma);
                }
                if (set_diag) {
                    size_t j;
                    for (j = 0; j < M; j++) { Ahicma[j * LDA + j] = diagVal; }
                }
                if (main_print_mat) {
                    printf("Ahicma\n");
                    printmat(Ahicma, M, M, LDA, MB, MB);
                }
                //LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', M, Ahicma, LDA);
                // AhicmaT: transpose of Ahicma
                PROGRESS("copy descAd into AhicmaT started");
                PASTE_TILE_TO_LAPACK(descAD, AhicmaT, check, double, LDA, M);

                {
                    size_t i, j;
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            Adense[j * LDA + i] = zero;
                        }
                    }
                }

                if (main_print_mat) {
                    printf("Ahicma-upperzero\n");
                    printmat(Ahicma, M, M, LDA, MB, MB);
                }
                PROGRESS("Transpose A started");
                LAPACKE_dge_trans(LAPACK_COL_MAJOR, M, M, Ahicma, LDA, AhicmaT, LDA);
                if (main_print_mat) {
                    printf("AhicmaT\n");
                    printmat(AhicmaT, M, M, LDA, MB, MB);
                }
                PROGRESS("TRMM started");
                cblas_dtrmm(CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, M, one, Ahicma, LDA,
                            AhicmaT, LDA);
                if (main_print_mat) {
                    printf("Ahicma*AhicmaT\n");
                    printmat(AhicmaT, M, M, LDA, MB, MB);
                }
                //double tmpnorm;normest(M, M, AhicmaT, &tmpnorm, swork);printf("tmpnorm:%e\n",tmpnorm);
                {
                    size_t i, j;
                    for (j = 0; j < M; j++) {
                        for (i = 0; i < j; i++) {
                            AhicmaT[j * LDA + i] = zero;
                        }
                    }
                }

                size_t nelm = M * M;
                if (main_print_mat)printf("nelm:%zu M:%d N:%d\n", nelm, M, N);
                PROGRESS("DAXPY started");
                cblas_daxpy(nelm, minusone, AhicmaT, 1, Adense, 1);
                if (main_print_mat) {
                    printf("Adense-(Ahicma*AhicmaT)\n");
                    printmat(Adense, M, M, LDA, MB, MB);
                }

                double normDenseAppDiff;
                PROGRESS("Norm of difference started");
                HCORE_dnormest(M, M, Adense, &normDenseAppDiff, swork);
                double accuracyDenseAppDiff = normDenseAppDiff / normA;
                //printf("normA:%.2e normDenseAppdiff:%.2e Accuracy: %.2e\n", normA, normDenseAppDiff,  accuracyDenseAppDiff);
                dparam[IPARAM_RES] = normDenseAppDiff;
                dparam[IPARAM_ANORM] = normA;
                dparam[IPARAM_XNORM] = normA;
                dparam[IPARAM_BNORM] = normAhicma;
#endif
            } else {
                PASTE_TILE_TO_LAPACK(descAD, Ahicma, check, double, LDA, M);
                PASTE_TILE_TO_LAPACK(descAD, AhicmaT, check, double, LDA, M);
            }
            PROGRESS("checking accuracy is finished");
        }
        filenumber = filenumber + (num_mpi_ranks);
        if (filenumber >= batches) {
            if (0)
                printf("\nI am node:%d and I am done filenumber:%d, batches:%f.\n", HICMA_My_Mpi_Rank(), filenumber,
                       batches);
            break;
        }

    } //batch for 
    STOP_TIMING();
    //free the starsh data structure 
    double maxcompute, mincompute, maxcompress, mincompress;
#ifdef HICMA_USE_MPI
    MPI_Allreduce( &total_compute_time, &maxcompute, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce( &total_compute_time, &mincompute, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

    MPI_Allreduce( &total_compress_time, &maxcompress, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce( &total_compress_time, &mincompress, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD);

   if(HICMA_My_Mpi_Rank()==0){
    printf("\nviruses: %d subviruses: %d, batches: %9.4f, batches_per_node: %d\n", iparam[IPARAM_NUMOBJ], iparam[IPARAM_NUMSUBOBJ], batches, batches_per_node);
    }
    if(HICMA_My_Mpi_Rank()==0){
    printf("\nMax Computation: %9.4f, Min Computation: %9.4f, Max Compression: %9.4f, Min Computation: %9.4f\n \n", maxcompute, mincompute, maxcompress, mincompress);
    }
#else
    printf("\nTLR Cholesky time: %9.4f, %9.4f \n\n", total_compute_time, total_compress_time);
#endif
    if (iparam[IPARAM_HICMA_STARSH_PROB] == PROBLEM_TYPE_3D_RBF_VIRUS) {
        STARSH_blrf *blrf = HICMA_get_starsh_format();
        STARSH_cluster *RC = blrf->row_cluster;
        void *RD = RC->data;
        starsh_mddata_free((STARSH_mddata *) RD);
    }

    PASTE_CODE_FREE_MATRIX(descAUV);
    PROGRESS("descAUV is freed");
    if (check == 0) { // If there is no check, then descAD and descDense are different. Refer to @KADIRLBL001
        PASTE_CODE_FREE_MATRIX(descAD);
        PROGRESS("descAD is freed");
    }
    PASTE_CODE_FREE_MATRIX(descArk);
    PROGRESS("descArk is freed");
    PASTE_CODE_FREE_MATRIX(descDense);
    if (solve)
        PASTE_CODE_FREE_MATRIX(descB);
    if (check_solve) {
        PASTE_CODE_FREE_MATRIX(descBcpy);
        PASTE_CODE_FREE_MATRIX(descmat);
    }
    PROGRESS("descDense is freed");
    PROGRESS("freed descs");
    return 0;
}
