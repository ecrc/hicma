/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

/**
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
/**
 * @file testing_dposv.c
 *
 *  HiCMA testing routine.
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2019-11-21
 * @brief This file shows how to factorize and solve using HiCMA. X/B matrix is in Tile Low Rank (TLR) format.
 **/
/*
 *
 * file testing_dposv.c
 *
 * HICMA testing routines
 * HICMA is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for HICMA 1.0.0
 * author Bilel Hadri, Hatem Ltaief
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * date 2010-11-15
 *
 */



#include "timing_auxiliary.h"
#include "timing.h"
#include "hicma_constants.h"
#include "hicma_struct.h"
#include "hicma_d.h"
#include "hicma.h"
#include "hicma_common.h"

#include "misc/auxcompute_d.h"
#include "misc/auxdescutil.h"

char datebuf[128];
time_t timer;
struct tm *tm_info;

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
int calc_rank_stat = 1;
int check = 0;


#include <time.h>
#include <sys/time.h>

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <hicma.h>
#include "testing_dauxiliary.h"

#if !defined(ARMPL)
enum blas_order_type {
    blas_rowmajor = 101,
    blas_colmajor = 102
};

enum blas_cmach_type {
    blas_base = 151,
    blas_t = 152,
    blas_rnd = 153,
    blas_ieee = 154,
    blas_emin = 155,
    blas_emax = 156,
    blas_eps = 157,
    blas_prec = 158,
    blas_underflow = 159,
    blas_overflow = 160,
    blas_sfmin = 161
};

enum blas_norm_type {
    blas_one_norm = 171,
    blas_real_one_norm = 172,
    blas_two_norm = 173,
    blas_frobenius_norm = 174,
    blas_inf_norm = 175,
    blas_real_inf_norm = 176,
    blas_max_norm = 177,
    blas_real_max_norm = 178
};

static void
BLAS_error(char *rname, int err, int val, int x) {
    fprintf(stderr, "%s %d %d %d\n", rname, err, val, x);
    abort();
}

#endif

static
void
BLAS_dge_norm(enum blas_order_type order, enum blas_norm_type norm,
              int m, int n, const double *a, int lda, double *res) {
    int i, j;
    float anorm, v;
    char rname[] = "BLAS_dge_norm";

    if (order != blas_colmajor) BLAS_error(rname, -1, order, 0);

    if (norm == blas_frobenius_norm) {
        anorm = 0.0f;
        for (j = n; j; --j) {
            for (i = m; i; --i) {
                v = a[0];
                anorm += v * v;
                a++;
            }
            a += lda - m;
        }
        anorm = sqrt(anorm);
    } else if (norm == blas_inf_norm) {
        anorm = 0.0f;
        for (i = 0; i < m; ++i) {
            v = 0.0f;
            for (j = 0; j < n; ++j) {
                v += cabs(a[i + j * lda]);
            }
            if (v > anorm)
                anorm = v;
        }
    } else {
        BLAS_error(rname, -2, norm, 0);
        return;
    }

    if (res) *res = anorm;
}

static
double
BLAS_dpow_di(double x, int n) {
    double rv = 1.0;

    if (n < 0) {
        n = -n;
        x = 1.0 / x;
    }

    for (; n; n >>= 1, x *= x) {
        if (n & 1)
            rv *= x;
    }

    return rv;
}

static
double
BLAS_dfpinfo(enum blas_cmach_type cmach) {
    double eps = 1.0, r = 1.0, o = 1.0, b = 2.0;
    int t = 53, l = 1024, m = -1021;
    char rname[] = "BLAS_dfpinfo";

    if ((sizeof eps) == sizeof(float)) {
        t = 24;
        l = 128;
        m = -125;
    } else {
        t = 53;
        l = 1024;
        m = -1021;
    }

    /* for (i = 0; i < t; ++i) eps *= half; */
    eps = BLAS_dpow_di(b, -t);
    /* for (i = 0; i >= m; --i) r *= half; */
    r = BLAS_dpow_di(b, m - 1);

    o -= eps;
    /* for (i = 0; i < l; ++i) o *= b; */
    o = (o * BLAS_dpow_di(b, l - 1)) * b;

    switch (cmach) {
        case blas_eps:
            return eps;
        case blas_sfmin:
            return r;
        default:
            BLAS_error(rname, -1, cmach, 0);
            break;
    }
    return 0.0;
}

int testing_dtrsmd(int argc, char **argv) {}//FIXME Use CMakeLists from Chameleon

static int check_factorization(int, double *, double *, int, int, double);

static int check_solution(int, int, double *, int, double *, double *, int, double);
/**
 * @name    Factorize and Solve
 * @brief   `timing_dposv()` shows how to use Cholesky factorization and then use triangular solve.
 * @ingroup testing
 */
/** This function shows matrix factorization and solve using TLR matrices.
 * 
 * Steps are as follows:
 */
int testing_dposv(int argc, char **argv) {

    //HICMA_set_print_index();
    //HICMA_set_print_index_end();
    //HICMA_unset_print_index_end();
    /*HICMA_set_use_fast_hcore_dgemm();*/

    int hres = 0;

    /*int nbnode = 0;*/
    /*HICMA_Comm_size( &nbnode );*/
    /*printf("Myrank:%d nbnode:%d\n", HICMA_My_Mpi_Rank(), nbnode);*/
    /*fflush( stdout );*/

    /* Check for number of arguments*/
    if (argc != 12) {
        USAGE("POSV", "N LDA NRHS LDB",
              "   - N    : the size of the matrix\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the RHS B\n"
              "   - NB   : Block size\n"
              "   - acc  : Fixed accuracy threshold\n"
              "   - rk   : Fixed rank threshold\n"
              "   - maxrank   : maxrank for U and V storage\n"
              "   - compmaxrank   : maxrank for starsh and buffers\n"
              "   - P   : num procs in row dim\n"
              "   - Q   : num procs in col dim\n"
              "   - check   : {0,1}\n"
        );
        return -1;
    }

    int N = atoi(argv[0]);
    int LDA = atoi(argv[1]);
    int NRHS = atoi(argv[2]);
    int LDB = atoi(argv[3]);
    int NB = atoi(argv[4]);
    double fixedacc = atof(argv[5]);
    int fixedrank = atoi(argv[6]);
    int maxrank = atoi(argv[7]);
    int comp_maxrank = atoi(argv[8]);
    int P = atoi(argv[9]);
    int Q = atoi(argv[10]);
    check = atoi(argv[11]);
    double eps;
    int uplo;
    int info_solution, info_factorization;
    int trans1, trans2;
    if (check) {
        printf("A: %d-by-%d LD:%d\n", N, N, LDA);
        printf("B: %d-by-%d LD:%d\n", N, NRHS, LDB);
        printf("NB: %d  fixedrank: %d  fixedacc:%.1e\n", NB, fixedrank, fixedacc);
        printf("MaxrankUV: %d  MaxrankBuffers: %d\n", maxrank, comp_maxrank);
    }

    double *A1 = NULL;
    double *A2 = NULL;
    double *B1 = NULL;
    double *B2 = NULL;

    if (check) {
        A1 = (double *) malloc(LDA * N * sizeof(double));
        A2 = (double *) malloc(LDA * N * sizeof(double));
        B1 = (double *) malloc(LDB * NRHS * sizeof(double));
        B2 = (double *) malloc(LDB * NRHS * sizeof(double));

        /* Check if unable to allocate memory */
        if ((!A1) || (!A2) || (!B1) || (!B2)) {
            printf("Out of Memory \n ");
            return -2;
        }
    }
    eps = BLAS_dfpinfo(blas_eps);

    uplo = HicmaUpper;
    uplo = HicmaLower; //to comply with current HICMA_dpotrf
    trans1 = uplo == HicmaUpper ? HicmaTrans : HicmaNoTrans;
    trans2 = uplo == HicmaUpper ? HicmaNoTrans : HicmaTrans;
    /*-------------------------------------------------------------
    *  TESTING DPOTRF + DPTRSM + DTRSM
    */

    /* Initialize A1 and A2 for Symmetric Positif Matrix */
    /*HICMA_dplgsy( N, HicmaUpperLower, N, A1, LDA, 51 );*/
    /*HICMA_dlacpy( HicmaUpperLower, N, N, A1, LDA, A2, LDA );*/

    /* Initialize B1 and B2 */
    /*HICMA_dplrnt( N, NRHS, B1, LDB, 371 );*/
    /*HICMA_dlacpy( HicmaUpperLower, N, NRHS, B1, LDB, B2, LDB );*/

    /*HICMA_dpotrf(uplo, N, A2, LDA);*/
    /*HICMA_dtrsm(HicmaLeft, uplo, trans1, HicmaNonUnit,*/
    /*N, NRHS, 1.0, A2, LDA, B2, LDB);*/
    /*HICMA_dtrsm(HicmaLeft, uplo, trans2, HicmaNonUnit,*/
    /*N, NRHS, 1.0, A2, LDA, B2, LDB);*/

    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 1;
    //eps = fixedacc/10;
    int NBxMaxrank = NB * maxrank;
    int NT = N / NB;
    int NTxMaxrank = NT * maxrank;
    int NBxNB = NB * NB;
    int NBxNRHS = NB * NRHS;
    int status;
    HICMA_problem_t hicma_problem;
    int probtype;
    char sym;
    double ddecay;
    /*BEGIN                AAAAAAAAAAAAAAAAAAAAAAAA                        */
    /** -# Allocate structures for TLR matrix \c A */
    HICMA_desc_t *descAUV = NULL;
    status = HICMA_Desc_Create(&descAUV, NULL, HicmaRealDouble, NB, maxrank * 2, NBxMaxrank * 2, LDA, NTxMaxrank * 2, 0,
                               0, N, NTxMaxrank * 2, P, Q);

    if (status != HICMA_SUCCESS) { return status; }
    HICMA_desc_t *descAD = NULL;

    status = HICMA_Desc_Create(&descAD, NULL, HicmaRealDouble, NB, NB, NBxNB, LDA, NB, 0, 0, N, NB, P, Q);
    if (status != HICMA_SUCCESS) { return status; }

    HICMA_desc_t *descArk = NULL;

    double *Ark_array = (double *) malloc(NT * NT * sizeof(double));

    memset(Ark_array, 0, NT * NT * sizeof(double));

    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'L', NT, NT, 0.0, NB, Ark_array, NT);

    status = HICMA_Desc_Create(&descArk, Ark_array, HicmaRealDouble, 1, 1, 1, NT, NT,
                               0, 0, NT, NT, P, Q);

    if (status != HICMA_SUCCESS) { return status; }

    HICMA_desc_t *descAdense = NULL;
    if (check) {
        status = HICMA_Desc_Create(&descAdense, NULL, HicmaRealDouble, NB, NB, NBxNB, LDA, N, 0, 0, N, N, P, Q);
        if (status != HICMA_SUCCESS) { return status; }
    }
    probtype = PROBLEM_TYPE_SS;
    sym = 'S';
    hicma_problem.ndim = 2;
    // Correlation length
    hicma_problem.beta = 0.1;
    // Smoothing parameter for Matern kernel
    hicma_problem.nu = 0.5;
    // Shift added to diagonal elements
    hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
    hicma_problem.noise = 5.e-4; //works for 640K but does not work for 10M
    hicma_problem.noise = 1.e-2; //

    // probtype = PROBLEM_TYPE_EDSIN;
    hicma_problem.wave_k = 20;
    hicma_problem.diag = N;

    generate_problem(probtype, sym, ddecay, N, NB, NT, NT, &hicma_problem);

    if (check) {
        /** -# Generate whole Dense Matrix \c descAdense for only checking purposes. */
        HICMA_dhagdm_Tile(HicmaLower, descAdense);
    }
    struct timeval tvalBefore, tvalAfter;
    gettimeofday(&tvalBefore, NULL);
    /** -# Generate diagonal Dense Matrix \c descAD  but put the tiles in a tall and skinny matrix. */
    HICMA_dhagdmdiag_Tile(HicmaUpperLower, descAD);
    gettimeofday(&tvalAfter, NULL);
    double tgad = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                  + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;
    gettimeofday(&tvalBefore, NULL);
    /** -# Generate Compressed Matrix \c descAUV. */
    HICMA_dhagcm_Tile(HicmaLower, descAUV, descArk, N, N, NB, NB, maxrank, fixedacc);
    gettimeofday(&tvalAfter, NULL);
    double tgauv = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                   + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;
    /** -# Check if there exists off diagonal dense blocks. */ //FIXME implement the function

    if (calc_rank_stat) { // print initial ranks of A
        PASTE_TILE_TO_LAPACK(descArk, Ark_initial, 1, double, NT, NT);
        if (HICMA_My_Mpi_Rank() == 0) {
            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaLower, Ark_initial, NT, NT, NT, &hicma_statrk_initial);
            printf("------------------------------A initial_ranks:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }
    double *Aorg = NULL;
    if (check) {
        PASTE_TILE_TO_LAPACK(descAdense, _Aorg, 1, double, LDA, N);
        Aorg = _Aorg;
        {  // fill upper part of Aorg so the rest of code taken from Chameleon works correctly.
            int i, j;
            for (j = 0; j < N; j++) {
                for (i = 0; i < j; i++) {
                    //Aorg[j*LDA+i] = 0.0;
                    Aorg[j * LDA + i] = Aorg[i * LDA + j];
                    //           Acham[j*LDA+i] = 0.0;
                }
            }
        }
    }
    /*END                AAAAAAAAAAAAAAAAAAAAAAAA                        */
    /*BEGIN                BBBBBBBBBBBBBBBBBBBBBBBB                        */
    /** -# Calculate number of tiles in column dimension of B via \f$ceil(NRHS*1.0/NB)\f$ */
    int NTB = ceil(NRHS * 1.0 / NB);
    int NTBxMaxrank = NTB * maxrank;
    HICMA_desc_t *descBUV = NULL;
    status = HICMA_Desc_Create(&descBUV, NULL, HicmaRealDouble, NB, maxrank * 2, NBxMaxrank * 2, LDA, NTBxMaxrank * 2,
                               0, 0, N, NTBxMaxrank * 2, P, Q);
    if (status != HICMA_SUCCESS) { return status; }
    HICMA_desc_t *descBrk = NULL;
    status = HICMA_Desc_Create(&descBrk, NULL, HicmaRealDouble, 1, 1, 1, NT, NTB, 0, 0, NT, NTB, P, Q);
    if (status != HICMA_SUCCESS) { return status; }
    HICMA_desc_t *descBdense = NULL;
    if (check) {
        status = HICMA_Desc_Create(&descBdense, NULL, HicmaRealDouble, NB, NB, NBxNB, LDA, NRHS, 0, 0, N, NRHS, P, Q);
        if (status != HICMA_SUCCESS) { return status; }
    }
    /*if randtlr*/
    //probtype = PROBLEM_TYPE_RND;
    sym = 'N';
    hicma_problem.noise = 0.0; //value added to diagonal
    //hicma_problem.wave_k = 20;
    hicma_problem.diag = hicma_problem.wave_k;

    ddecay = 1e-2;
    double fixedacc_B = 1e-12;//fixedacc;
    fixedacc_B = fixedacc;
    generate_problem(probtype, sym, ddecay, N, NB, NT, NT, &hicma_problem);

    if (check) {
        /** -# Generate whole Dense Matrix \c descBdense for only checking purposes. */
        HICMA_dhagdm_Tile(HicmaUpperLower, descBdense);
    }
    /** -# Generate Compressed Matrix \c descBUV. */
    gettimeofday(&tvalBefore, NULL);
    HICMA_dhagcm_Tile(HicmaUpperLower, descBUV, descBrk, N, NRHS, NB, NB, maxrank, fixedacc_B);
    gettimeofday(&tvalAfter, NULL);
    double tgbuv = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                   + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;
    //check_same(descBrk, descBrk2, 'A', 'A');
    /** -# Check if there exists off diagonal dense blocks. */ //FIXME implement the function

    if (calc_rank_stat) { // print initial ranks
        PASTE_TILE_TO_LAPACK(descBrk, Brk_initial, 1, double, NT, NTB);
        if (HICMA_My_Mpi_Rank() == 0) {
            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaLower, Brk_initial, NT, NTB, NT, &hicma_statrk_initial);
            printf("------------------------------B initial_ranks:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }
    double *Borg = NULL;
    if (check) {
        PASTE_TILE_TO_LAPACK(descBdense, _Borg, 1, double, LDA, NRHS);
        Borg = _Borg;
        /*END                BBBBBBBBBBBBBBBBBBBBBBBB                        */
        // Save A
//        HICMA_dlacpy(HicmaUpperLower, N, N, Aorg, LDA, A1, LDA);

        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'X', N, N, Aorg, LDA, A1, LDA);

//        HICMA_dlacpy(HicmaUpperLower, N, N, Aorg, LDA, A2, LDA);
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'X', N, N, Aorg, LDA, A2, LDA);

        if (main_print_mat) {
            printf("Aorg\n");
            printmat_format(Aorg, N, N, LDA, 1000, 1000, 0);
        }
        // Save B
//        HICMA_dlacpy(HicmaUpperLower, N, NRHS, Borg, LDA, B1, LDA);
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'X', N, NRHS, Borg, LDA, B1, LDA);
//        HICMA_dlacpy(HicmaUpperLower, N, NRHS, Borg, LDA, B2, LDA);
        LAPACKE_dlacpy(LAPACK_COL_MAJOR, 'X', N, NRHS, Borg, LDA, B2, LDA);
        if (main_print_mat) {
            printf("Borg\n");
            printmat_format(Borg, N, NRHS, LDB, 1000, 1000, 0);
        }

//        HICMA_dpotrf(uplo, N, A2, LDA);
        char uplo_c;
        CBLAS_UPLO uplo_cblas;
        if (uplo == HicmaUpper) {
            uplo_c = 'U';
            uplo_cblas = CblasUpper;
        } else {
            uplo_c = 'L';
            uplo_cblas = CblasLower;
        }

        LAPACKE_dpotrf(LAPACK_COL_MAJOR, uplo_c, N, A2, LDA);

        CBLAS_TRANSPOSE trans1_cblas;
        if (trans1 == HicmaTrans) {
            trans1_cblas = CblasTrans;
        } else {
            trans1_cblas = CblasNoTrans;
        }

        cblas_dtrsm(CblasColMajor, CblasLeft, uplo_cblas, trans1_cblas, CblasNonUnit, N, NRHS, 1.0, A2, LDA, B2, LDB);

//        HICMA_dtrsm(HicmaLeft, uplo, trans1, HicmaNonUnit,
//                    N, NRHS, 1.0, A2, LDA, B2, LDB);

        CBLAS_TRANSPOSE trans2_cblas;
        if (trans2 == HicmaTrans) {
            trans2_cblas = CblasTrans;
        } else {
            trans2_cblas = CblasNoTrans;
        }
        cblas_dtrsm(CblasColMajor, CblasLeft, uplo_cblas, trans2_cblas, CblasNonUnit, N, NRHS, 1.0, A2, LDA, B2, LDB);

//        HICMA_dtrsm(HicmaLeft, uplo, trans2, HicmaNonUnit,
//                    N, NRHS, 1.0, A2, LDA, B2, LDB);
    }
    /** -# Perform Cholesky factorization. */
    gettimeofday(&tvalBefore, NULL);
    HICMA_dpotrf_Tile(HicmaLower, descAUV, descAD, descArk, fixedrank, comp_maxrank, fixedacc);
    gettimeofday(&tvalAfter, NULL);
    double tpotrf = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                    + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;
    double *Ahicma = NULL;
    if (check) {
        /** -# Uncompress low rank tiles. */
        HICMA_duncompress(HicmaLower, descAUV, descAdense, descArk);
        /** -# Scatter dense tiles on to diagonal of dense matrix */
        HICMA_ddiag_vec2mat(descAD, descAdense);
        PASTE_TILE_TO_LAPACK(descAdense, _Ahicma, 1, double, LDA, N);
        Ahicma = _Ahicma;
    }
    if (calc_rank_stat) { // print final ranks of L
        PASTE_TILE_TO_LAPACK(descArk, Ark_initial, 1, double, NT, NT);
        if (HICMA_My_Mpi_Rank() == 0) {
            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaLower, Ark_initial, NT, NT, NT, &hicma_statrk_initial);
            printf("------------------------------A final_ranks after HICMA_potrf():");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }

    if (check) {
        //Check Chameleon. A1 is original matrix, A2 is L coming from HICMA_dpotrf
        //info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
        //Check Hicma. Aorg is original matrix, Ahicma is L coming from HICMA_dpotrf
        info_factorization = check_factorization(N, Aorg, Ahicma, LDA, uplo, fixedacc);
        //return info_factorization; 

        //check if uncompress works well for B
        /*HICMA_duncompress(HicmaUpperLower, descBUV, descBdense2, descBrk); */
        /*PASTE_TILE_TO_LAPACK( descBdense2, Bhicma2, 1, double, LDA, N );*/
        /*HICMA_dtrsm(HicmaLeft, uplo, trans1, HicmaNonUnit,*/
        /*N, NRHS, 1.0, A2, LDA, Bhicma2, LDB);*/
        /*HICMA_dtrsm(HicmaLeft, uplo, trans2, HicmaNonUnit,*/
        /*N, NRHS, 1.0, A2, LDA, Bhicma2, LDB);*/
        /*info_solution = check_solution(N, NRHS, Aorg, LDA, Borg, Bhicma2, LDB, eps);*/
        /*return 0;*/

        printf("\n");
        printf("------ TESTS FOR CHAMELEON DPOTRF + DTRSM + DTRSM  ROUTINE -------  \n");
        printf("            Size of the Matrix %d by %d\n", N, N);
        printf("\n");
        printf(" The matrix A is randomly generated for each test.\n");
        printf("============\n");
        printf(" The relative machine precision (eps) is to be %e \n", eps);
        printf(" Computational tests pass if scaled residuals are less than 60.\n");
    }

    gettimeofday(&tvalBefore, NULL);

    HICMA_dtrsm_Tile(HicmaLeft, uplo, trans1, HicmaNonUnit, 1.0,
                     descAUV,
                     descAD,
                     descArk,
                     descBUV,
                     descBrk,
                     fixedrank, comp_maxrank, fixedacc);
    gettimeofday(&tvalAfter, NULL);
    double ttrsm1 = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                    + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;

    if (calc_rank_stat) { // print ranks for B after 1st trsm
        PASTE_TILE_TO_LAPACK(descBrk, Brk_initial, 1, double, NT, NTB);
        if (HICMA_My_Mpi_Rank() == 0) {
            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaUpperLower, Brk_initial, NT, NTB, NT, &hicma_statrk_initial);
            printf("------------------------------B ranks after 1st trsm:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }
    gettimeofday(&tvalBefore, NULL);

    HICMA_dtrsm_Tile(HicmaLeft, uplo, trans2, HicmaNonUnit, 1.0,
                     descAUV,
                     descAD,
                     descArk,
                     descBUV,
                     descBrk,
                     fixedrank, comp_maxrank, fixedacc);
    gettimeofday(&tvalAfter, NULL);
    double ttrsm2 = (tvalAfter.tv_sec - tvalBefore.tv_sec)
                    + (tvalAfter.tv_usec - tvalBefore.tv_usec) / 1000000.0;

    if (calc_rank_stat) { // print initial ranks
        PASTE_TILE_TO_LAPACK(descBrk, Brk_initial, 1, double, NT, NTB);
        if (HICMA_My_Mpi_Rank() == 0) {
            HICMA_stat_t hicma_statrk_initial;
            dget_stat(HicmaUpperLower, Brk_initial, NT, 1, NT, &hicma_statrk_initial);
            printf("------------------------------B ranks after 2nd trsm:");
            dprint_stat(hicma_statrk_initial);
            fflush(stderr);
            fflush(stdout);
        }
    }

    if (HICMA_My_Mpi_Rank() == 0) {
        printf("%d %d %d %d %d %d %d %d %d %.1e %d %d ", P * Q, P, Q, N, NRHS, LDA, LDB, NB, fixedrank, fixedacc,
               maxrank, comp_maxrank);
        printf("%g %g %g %g %g %g\n", tgad, tgauv, tgbuv, tpotrf, ttrsm1, ttrsm2);
        fflush(stderr);
        fflush(stdout);
    }

    if (check) {
        //printf("%s %d %s EXITING %e\n", __FILE__, __LINE__, __func__, fixedacc);exit(0);
        HICMA_duncompress_custom_size(HicmaUpperLower, descBUV, descBdense, descBrk, N, NRHS, NB, NB);
        PASTE_TILE_TO_LAPACK(descBdense, Bhicma, 1, double, LDA, NRHS);
        info_solution = check_solution(N, NRHS, Aorg, LDA, Borg, Bhicma, LDB, fixedacc);
        /* Check the factorization and the solution */
        //info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
        info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

        if ((info_solution == 0) & (info_factorization == 0)) {
            printf("***************************************************\n");
            printf(" ---- TESTING DPOTRF + DTRSM + DTRSM ..... PASSED !\n");
            printf("***************************************************\n");
        } else {
            printf("***************************************************\n");
            printf(" - TESTING DPOTRF + DTRSM + DTRSM ... FAILED !\n");
            printf("***************************************************\n");
        }
        free(A1);
        free(A2);
        free(B1);
        free(B2);
    }
    if (check) {
        PASTE_CODE_FREE_MATRIX(descAdense);
        PASTE_CODE_FREE_MATRIX(descBdense);
    }

    PASTE_CODE_FREE_MATRIX(descAD);
    PASTE_CODE_FREE_MATRIX(descAUV);
    PASTE_CODE_FREE_MATRIX(descArk);

    PASTE_CODE_FREE_MATRIX(descBUV);
    PASTE_CODE_FREE_MATRIX(descBrk);

    return hres;
}


/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
static int check_factorization(int N, double *A1, double *A2, int LDA, int uplo, double eps) {
    double Anorm, Rnorm;
    double alpha;
    int info_factorization;
    int i, j;

    double *Residual = (double *) malloc(N * N * sizeof(double));
    double *L1 = (double *) malloc(N * N * sizeof(double));
    double *L2 = (double *) malloc(N * N * sizeof(double));
    double *work = (double *) malloc(N * sizeof(double));

    memset((void *) L1, 0, N * N * sizeof(double));
    memset((void *) L2, 0, N * N * sizeof(double));

    alpha = 1.0;

    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, ' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == HicmaUpper) {
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'u', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'u', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    } else {
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'l', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR, 'l', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            Residual[j * N + i] = L2[j * N + i] - Residual[j * N + i];

    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, N, Residual, N, &Rnorm);
    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, N, A1, LDA, &Anorm);

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e eps:%.1e\n", Rnorm / (Anorm * N * eps), eps);

    if (isnan(Rnorm / (Anorm * N * eps)) || isinf(Rnorm / (Anorm * N * eps)) || (Rnorm / (Anorm * N * eps) > 60.0)) {
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    } else {
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual);
    free(L1);
    free(L2);
    free(work);

    return info_factorization;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

static int check_solution(int N, int NRHS, double *A1, int LDA, double *B1, double *B2, int LDB, double eps) {
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double alpha, beta;
    double *work = (double *) malloc(N * sizeof(double));

    alpha = 1.0;
    beta = -1.0;

    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, NRHS, B2, LDB, &Xnorm);
    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, N, A1, LDA, &Anorm);
    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, NRHS, B1, LDB, &Bnorm);

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, (alpha), A1, LDA, B2, LDB, (beta), B1, LDB);
    BLAS_dge_norm(blas_colmajor, blas_inf_norm, N, NRHS, B1, LDB, &Rnorm);

    if (getenv("HiCMA_TESTING_VERBOSE")) {
        printf("||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e eps:%.2e\n", Anorm, Xnorm, Bnorm, Rnorm, eps);
    }
    result = Rnorm / ((Anorm * Xnorm + Bnorm) * N * eps);
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0)) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
    } else {
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
