/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/
/**
 * @file testing_ztrsmd.c
 *
 *  HiCMA testing routine.
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @brief This file shows how to factorize and solve using HiCMA. X/B matrix in AX=B  is dense.
 **/
/*
 *
 * file testing_zposv.c
 *
 * MORSE testing routines
 * MORSE is a software package provided by Univ. of Tennessee,
 * Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * version 2.5.0
 * comment This file has been automatically generated
 *         from Plasma 2.5.0 for MORSE 1.0.0
 * author Bilel Hadri, Hatem Ltaief
 * author Mathieu Faverge
 * author Emmanuel Agullo
 * author Cedric Castagnede
 * date 2010-11-15
 * precisions normal z -> c d s
 *
 */



#include "timing.h"
#include "hicma_constants.h"
#include "hicma_struct.h"
#include "hicma_z.h"
#include "hicma.h"
#include "hicma_common.h"

#include <assert.h>
#include "auxcompute_z.h"
#include "auxdescutil.h"

// zgytlr uses starsh in MPI mode.
STARSH_blrf *mpiF;

int print_progress = 1;   // Print progress about the execution
char datebuf[128];
time_t timer;
struct tm* tm_info;

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


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include <morse.h>
#include "testing_zauxiliary.h"

enum blas_order_type {
            blas_rowmajor = 101,
            blas_colmajor = 102 };

enum blas_cmach_type {
            blas_base      = 151,
            blas_t         = 152,
            blas_rnd       = 153,
            blas_ieee      = 154,
            blas_emin      = 155,
            blas_emax      = 156,
            blas_eps       = 157,
            blas_prec      = 158,
            blas_underflow = 159,
            blas_overflow  = 160,
            blas_sfmin     = 161};

enum blas_norm_type {
            blas_one_norm       = 171,
            blas_real_one_norm  = 172,
            blas_two_norm       = 173,
            blas_frobenius_norm = 174,
            blas_inf_norm       = 175,
            blas_real_inf_norm  = 176,
            blas_max_norm       = 177,
            blas_real_max_norm  = 178 };

static void
BLAS_error(char *rname, int err, int val, int x) {
  fprintf( stderr, "%s %d %d %d\n", rname, err, val, x );
  abort();
}

static
void
BLAS_dge_norm(enum blas_order_type order, enum blas_norm_type norm,
  int m, int n, const double *a, int lda, double *res) {
  int i, j; float anorm, v;
  char rname[] = "BLAS_dge_norm";

  if (order != blas_colmajor) BLAS_error( rname, -1, order, 0 );

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
    anorm = sqrt( anorm );
  } else if (norm == blas_inf_norm) {
    anorm = 0.0f;
    for (i = 0; i < m; ++i) {
      v = 0.0f;
      for (j = 0; j < n; ++j) {
        v += cabs( a[i + j * lda] );
      }
      if (v > anorm)
        anorm = v;
    }
  } else {
    BLAS_error( rname, -2, norm, 0 );
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
  eps = BLAS_dpow_di( b, -t );
  /* for (i = 0; i >= m; --i) r *= half; */
  r = BLAS_dpow_di( b, m-1 );

  o -= eps;
  /* for (i = 0; i < l; ++i) o *= b; */
  o = (o * BLAS_dpow_di( b, l-1 )) * b;

  switch (cmach) {
    case blas_eps: return eps;
    case blas_sfmin: return r;
    default:
      BLAS_error( rname, -1, cmach, 0 );
      break;
  }
  return 0.0;
}

int testing_dposv(int argc, char** argv){}//FIXME Use CMakeLists from Chameleon

static int check_factorization(int, double*, double*, int, int , double);
static int check_solution(int, int, double*, int, double*, double*, int, double);
/**
 * @name    Factorize and Solve
 * @brief   `timing_dtrsmd()` shows how to use Cholesky factorization and then use triangular solve.
 * @ingroup testing
 */
/** This function shows matrix factorization and solve using TLR and dense matrices.
 * A is in TLR format. X/B is in dense format.
 * 
 * Steps are as follows:
 */
int testing_dtrsmd(int argc, char **argv)
{
    HICMA_init();
    /*HICMA_set_print_index();*/
    /*HICMA_set_print_mat();*/
    //HICMA_set_print_index_end();
    //HICMA_unset_print_index_end();
    /*HICMA_set_use_fast_hcore_zgemm();*/
    
    int hres = 0;

    /* Check for number of arguments*/
    if (argc != 7){
        USAGE("TRSMD", "N LDA NRHS LDB",
              "   - N    : the size of the matrix\n"
              "   - LDA  : leading dimension of the matrix A\n"
              "   - NRHS : number of RHS\n"
              "   - LDB  : leading dimension of the RHS B\n"
              "   - NB   : Block size\n"
              "   - acc  : Fixed accuracy threshold\n"
              "   - rk   : Fixed rank threshold\n"
              );
        return -1;
    }

    int N     = atoi(argv[0]);
    int LDA   = atoi(argv[1]);
    int NRHS  = atoi(argv[2]);
    int LDB   = atoi(argv[3]);
    int NB    = atoi(argv[4]);
    double fixedacc    = atof(argv[5]);
    int fixedrank    = atoi(argv[6]);
    double eps;
    int uplo;
    int info_solution, info_factorization;
    int trans1, trans2;
    printf("A: %d-by-%d LD:%d\n", N, N, LDA);
    printf("B: %d-by-%d LD:%d\n", N, NRHS, LDB);
    printf("NB: %d  fixedrank: %d  fixedacc:%.1e\n", NB ,fixedrank, fixedacc);

    double *A1   = (double *)malloc(LDA*N*sizeof(double));
    double *A2   = (double *)malloc(LDA*N*sizeof(double));
    double *B1   = (double *)malloc(LDB*NRHS*sizeof(double));
    double *B2   = (double *)malloc(LDB*NRHS*sizeof(double));

    /* Check if unable to allocate memory */
    if ((!A1)||(!A2)||(!B1)||(!B2)){
        printf("Out of Memory \n ");
        return -2;
    }

    eps = BLAS_dfpinfo( blas_eps );

    uplo = MorseUpper;
    uplo = MorseLower; //to comply with current HICMA_dpotrf
    trans1 = uplo == MorseUpper ? MorseTrans : MorseNoTrans;
    trans2 = uplo == MorseUpper ? MorseNoTrans : MorseTrans;
    /*-------------------------------------------------------------
    *  TESTING ZPOTRF + ZPTRSM + ZTRSM
    */

    /* Initialize A1 and A2 for Symmetric Positif Matrix */
    /*MORSE_dplgsy( N, MorseUpperLower, N, A1, LDA, 51 );*/
    /*MORSE_dlacpy( MorseUpperLower, N, N, A1, LDA, A2, LDA );*/

    /* Initialize B1 and B2 */
    /*MORSE_dplrnt( N, NRHS, B1, LDB, 371 );*/
    /*MORSE_dlacpy( MorseUpperLower, N, NRHS, B1, LDB, B2, LDB );*/

    /*MORSE_dpotrf(uplo, N, A2, LDA);*/
    /*MORSE_dtrsm(MorseLeft, uplo, trans1, MorseNonUnit,*/
                 /*N, NRHS, 1.0, A2, LDA, B2, LDB);*/
    /*MORSE_dtrsm(MorseLeft, uplo, trans2, MorseNonUnit,*/
                 /*N, NRHS, 1.0, A2, LDA, B2, LDB);*/
    
	// DO NOT enforce compression of diagonal tiles, required for HICMA_zgytlr()
	int compress_diag = 0;
    global_check = 1;
    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 1;
    int P = 1;
    int Q = 1;
    int maxrank = NB/2;
    int comp_maxrank = NB;
    //eps = fixedacc/10;
    int NBxMaxrank = NB * maxrank;
    int NT = 0;
    if (N % NB == 0)
        NT = N / NB;
    else
        NT = N / NB + 1;
    int NTxMaxrank = NT * maxrank;
    int NBxNB = NB * NB;
    int NBxNRHS = NB * NRHS;
    int NTxNB = NT * NB;
    int status;
    HICMA_problem_t hicma_problem;
    int probtype;
    char sym;
    double ddecay;
    /*BEGIN                AAAAAAAAAAAAAAAAAAAAAAAA                        */
    /** -# Allocate structures for TLR matrix \c A */
    MORSE_desc_t *descAUV = NULL;
    status = MORSE_Desc_Create(&descAUV, NULL, MorseRealDouble, NB, maxrank*2, NBxMaxrank*2, NTxNB,  NTxMaxrank*2, 0, 0, NTxNB, NTxMaxrank*2, P, Q);
    printf("AUV: m:%d n:%d mb:%d nb:%d mt:%d nt:%d\n", 
            descAUV->m,
            descAUV->n,
            descAUV->mb,
            descAUV->nb,
            descAUV->mt,
            descAUV->nt
            );
    if(status != MORSE_SUCCESS){ return status; } 
    MORSE_desc_t *descAD = NULL;
    status = MORSE_Desc_Create(&descAD, NULL, MorseRealDouble, NB, NB, NBxNB, LDA,  NB, 0, 0, N, NB, P, Q); 
    if(status != MORSE_SUCCESS){ return status; } 
    
    MORSE_desc_t *descArk = NULL;
    status = MORSE_Desc_Create(&descArk, NULL, MorseRealDouble, 1, 1, 1, NT, NT, 0, 0, NT, NT, P, Q);
    if(status != MORSE_SUCCESS){ return status; } 
    /** -# Allocate rank matrix for \c A. Set off-diagonal elements to 0 and diagonal elements to NB since only diagonals are dense. */
    MORSE_dlaset_Tile(MorseLower, 0.0, NB, descArk);

    MORSE_desc_t *descAdense = NULL;
    status = MORSE_Desc_Create(&descAdense, NULL, MorseRealDouble, NB, NB, NBxNB, LDA, N, 0, 0, N, N, P, Q);
    if(status != MORSE_SUCCESS){ return status; } 
    // For checking against chameleon 
    MORSE_desc_t *descAdense2 = NULL;
    status = MORSE_Desc_Create(&descAdense2, NULL, MorseRealDouble, NB, NB, NBxNB, LDA, N, 0, 0, N, N, P, Q);
    if(status != MORSE_SUCCESS){ return status; } 
    probtype = HICMA_STARSH_PROB_SS;
    sym = 'S';
    hicma_problem.ndim = 2;
    // Correlation length
    hicma_problem.beta  = 0.1;
    // Smoothing parameter for Matern kernel
    hicma_problem.nu    = 0.5;
    // Shift added to diagonal elements
    hicma_problem.noise = 1.e-4; //not enough for matrices larger than 600K
    hicma_problem.noise = 5.e-4; //works for 640K but does not work for 10M
    hicma_problem.noise = 1.e-2; //

    HICMA_zgenerate_problem(probtype, sym, ddecay, N, NB, NT, NT, &hicma_problem);

    printf("Generate dense A\n");
    /** -# Generate whole Dense Matrix \c descAdense for only checking purposes. */
    HICMA_zhagdm_Tile(MorseLower, descAdense);
    printf("Generate diagonals of A\n");
    /** -# Generate diagonal Dense Matrix \c descAD  but put the tiles in a tall and skinny matrix. */
    HICMA_zhagdmdiag_Tile(MorseUpperLower, descAD);
    printf("Generate low rank off-diagonals of A\n");
    /** -# Generate Compressed Matrix \c descAUV. */
    status = HICMA_zhagcm_Tile(MorseLower, descAUV, descArk, N, N, NB, NB, maxrank, fixedacc);
    if(status != MORSE_SUCCESS){ return status; } 
    /** -# Check if there exists off diagonal dense blocks. */ //FIXME implement the function

    if(calc_rank_stat){ // print initial ranks of A
        PASTE_TILE_TO_LAPACK( descArk, Ark_initial, 1, double, NT, NT );
        HICMA_stat_t hicma_statrk_initial;
        zget_stat(MorseLower, Ark_initial, NT, NT, NT,  &hicma_statrk_initial);
        if(MORSE_My_Mpi_Rank()==0){
            printf("------------------------------A initial_ranks:");
            zprint_stat(hicma_statrk_initial);
        }
        fflush(stderr);
        fflush(stdout);
    }
    PASTE_TILE_TO_LAPACK( descAdense, Aorg, 1, double, LDA, N );
    {  // fill upper part of Aorg so the rest of code taken from Chameleon works correctly.
        int i, j;
        for(j = 0; j < N; j++){
            for(i = 0; i < j; i++){
                //Aorg[j*LDA+i] = 0.0;
                Aorg[j*LDA+i] = Aorg[i*LDA+j];
     //           Acham[j*LDA+i] = 0.0;
            }
        }
    }
    /*END                AAAAAAAAAAAAAAAAAAAAAAAA                        */
    /*BEGIN                BBBBBBBBBBBBBBBBBBBBBBBB                        */
    /** -# Calculate number of tiles in column dimension of B via \f$ceil(NRHS*1.0/NB)\f$ */
    int NTB = ceil(NRHS*1.0/NB);
    MORSE_desc_t *descBdense = NULL;
    status = MORSE_Desc_Create(&descBdense, NULL, MorseRealDouble, NB, NB, NBxNB, LDA, NRHS, 0, 0, N, NRHS, P, Q);
    if(status != MORSE_SUCCESS){ return status; } 

    // For checking against chameleon 
    MORSE_desc_t *descBdense2 = NULL;
    status = MORSE_Desc_Create(&descBdense2, NULL, MorseRealDouble, NB, NB, NBxNB, LDA, NRHS, 0, 0, N, NRHS, P, Q);
    if(status != MORSE_SUCCESS){ return status; } 
    /*if randtlr*/
    probtype = HICMA_STARSH_PROB_RND;
    sym = 'N';
    hicma_problem.noise    = 0.0; //value added to diagonal

    ddecay = 1e-2;
    //double fixedacc_B = 1e-3;
    double fixedacc_B = fixedacc;
    //printf("Generate problem for B\n");
    //HICMA_zgenerate_problem(probtype, sym, ddecay, N, NB, NT, NT, &hicma_problem);
    //mpiF = hicma_problem.starsh_format; // This is assignment will be hidden from user in release

    printf("Generate dense B\n");
    MORSE_dplrnt_Tile(descBdense, 371 );
    MORSE_dlacpy_Tile( MorseUpperLower, descBdense, descBdense2 );
    //MORSE_dplrnt_Tile(descBdense2, 371 );
    /** -# Generate whole Dense Matrix \c descBdense for only checking purposes. */
    //HICMA_zhagdm_Tile(MorseUpperLower, descBdense);
    PASTE_TILE_TO_LAPACK( descBdense, Borg, 1, double, LDA, NRHS ); 
    /*END                BBBBBBBBBBBBBBBBBBBBBBBB                        */
    // Save A
    MORSE_dlacpy( MorseUpperLower, N, N, Aorg, LDA, A1, LDA );
    MORSE_dlacpy( MorseUpperLower, N, N, Aorg, LDA, A2, LDA );
    if(main_print_mat ){printf("Aorg\n");printmat_format(Aorg,N,N,LDA,1000, 1000, 0);}
    // Save B
    MORSE_dlacpy( MorseUpperLower, N, NRHS, Borg, LDA, B1, LDA );
    MORSE_dlacpy( MorseUpperLower, N, NRHS, Borg, LDA, B2, LDA );
    if(main_print_mat ){printf("Borg\n");printmat_format(Borg,N,NRHS,LDB,1000, 1000, 0);}
    
    
    MORSE_dpotrf(uplo, N, A2, LDA);
    MORSE_dtrsm(MorseLeft, uplo, trans1, MorseNonUnit,
                 N, NRHS, 1.0, A2, LDA, B2, LDB);
    if(0)cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, 1.0, A1, LDA, B2, LDB, 0.0, B1, LDB);
    MORSE_dtrsm(MorseLeft, uplo, trans2, MorseNonUnit,
                 N, NRHS, 1.0, A2, LDA, B2, LDB);
    
    printf("Potrf on A\n");
    /** -# Perform Cholesky factorization. */
    HICMA_zpotrf_Tile(MorseLower, descAUV, descAD, descArk, fixedrank, comp_maxrank, fixedacc );
    /** -# Uncompress low rank tiles. */
    HICMA_zuncompress(MorseLower, descAUV, descAdense, descArk);
    /** -# Scatter dense tiles on to diagonal of dense matrix */
    HICMA_zdiag_vec2mat(descAD, descAdense);
    PASTE_TILE_TO_LAPACK( descAdense, Ahicma, 1, double, LDA, N );
    if(calc_rank_stat){ // print final ranks of L
        PASTE_TILE_TO_LAPACK( descArk, Ark_initial, 1, double, NT, NT ); 
        HICMA_stat_t hicma_statrk_initial;
        zget_stat(MorseLower, Ark_initial, NT, NT, NT,  &hicma_statrk_initial);
        if(MORSE_My_Mpi_Rank()==0){
            printf("------------------------------A final_ranks after HICMA_potrf():");
            zprint_stat(hicma_statrk_initial);
        }
        fflush(stderr);
        fflush(stdout);
    }

    //Check Chameleon. A1 is original matrix, A2 is L coming from MORSE_dpotrf
    //info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
    //Check Hicma. Aorg is original matrix, Ahicma is L coming from HICMA_dpotrf
    info_factorization = check_factorization( N, Aorg, Ahicma, LDA, uplo, fixedacc);
    //return info_factorization; 

    //check if uncompress works well for B
    /*HICMA_zuncompress(MorseUpperLower, descBUV, descBdense2, descBrk); */
    /*PASTE_TILE_TO_LAPACK( descBdense2, Bhicma2, 1, double, LDA, N );*/
    /*MORSE_dtrsm(MorseLeft, uplo, trans1, MorseNonUnit,*/
                 /*N, NRHS, 1.0, A2, LDA, Bhicma2, LDB);*/
    /*MORSE_dtrsm(MorseLeft, uplo, trans2, MorseNonUnit,*/
                 /*N, NRHS, 1.0, A2, LDA, Bhicma2, LDB);*/
    /*info_solution = check_solution(N, NRHS, Aorg, LDA, Borg, Bhicma2, LDB, eps);*/
    /*return 0;*/

    printf("\n");
    printf("------ TESTS FOR CHAMELEON ZPOTRF + ZTRSM + ZTRSM  ROUTINE -------  \n");
    printf("            Size of the Matrix %d by %d\n", N, N);
    printf("\n");
    printf(" The matrix A is randomly generated for each test.\n");
    printf("============\n");
    printf(" The relative machine precision (eps) is to be %e \n", eps);
    printf(" Computational tests pass if scaled residuals are less than 60.\n");

    PASTE_TILE_TO_LAPACK( descBdense, arrayB_1, 1, double, LDA, NRHS ); 
    PASTE_TILE_TO_LAPACK( descBdense2, arrayB2_1, 1, double, LDA, NRHS ); 
    check_same_array(arrayB_1, arrayB2_1, N, __LINE__, __FILE__);

    printf("First trsm\n");
    if(1)MORSE_dtrsm_Tile(MorseLeft, uplo, trans1, MorseNonUnit, 1.0, 
            descAdense, descBdense2);
    if(1)HICMA_ztrsmd_Tile(MorseLeft, uplo, trans1, MorseNonUnit, 1.0, 
            descAUV, 
            descAD, 
            descArk,
            descBdense, comp_maxrank);
    PASTE_TILE_TO_LAPACK( descBdense, arrayB_2, 1, double, LDA, NRHS ); 
    PASTE_TILE_TO_LAPACK( descBdense2, arrayB2_2, 1, double, LDA, NRHS ); 
    check_same_array(arrayB_2, arrayB2_2, N, __LINE__, __FILE__);

    printf("Second trsm\n");
    if(1)MORSE_dtrsm_Tile(MorseLeft, uplo, trans2, MorseNonUnit, 1.0, 
            descAdense, descBdense2);
    if(1)HICMA_ztrsmd_Tile(MorseLeft, uplo, trans2, MorseNonUnit, 1.0, 
            descAUV, 
            descAD, 
            descArk,
            descBdense, comp_maxrank);
    PASTE_TILE_TO_LAPACK( descBdense, arrayB_3, 1, double, LDA, NRHS ); 
    PASTE_TILE_TO_LAPACK( descBdense2, arrayB2_3, 1, double, LDA, NRHS ); 
    check_same_array(arrayB_3, arrayB2_3, N, __LINE__, __FILE__);

    //printf("%s %d %s EXITING %e\n", __FILE__, __LINE__, __func__, fixedacc);exit(0);
    PASTE_TILE_TO_LAPACK( descBdense, Bhicma, 1, double, LDA, NRHS );
    int info_solution_hicma = check_solution(N, NRHS, Aorg, LDA, Borg, Bhicma, LDB, fixedacc);
    //info_solution = check_solution(N, NRHS, Aorg, LDA, Borg, Bhicma, LDB, eps);
    /* Check the factorization and the solution */
    //info_factorization = check_factorization( N, A1, A2, LDA, uplo, eps);
    info_solution = check_solution(N, NRHS, A1, LDA, B1, B2, LDB, eps);

    if ((info_solution_hicma == 0)&(info_solution == 0)&(info_factorization == 0)){
        printf("***************************************************\n");
        printf(" ---- TESTING ZPOTRF + ZTRSM + ZTRSM ..... PASSED !\n");
        printf("***************************************************\n");
    }
    else{
        printf("***************************************************\n");
        printf(" - TESTING ZPOTRF + ZTRSM + ZTRSM ... FAILED !\n");
        printf("***************************************************\n");
    }

    free(A1); free(A2); free(B1); free(B2);
    PASTE_CODE_FREE_MATRIX( descAdense );
    PASTE_CODE_FREE_MATRIX( descBdense );
    PASTE_CODE_FREE_MATRIX( descAdense2 );
    PASTE_CODE_FREE_MATRIX( descBdense2 );
    PASTE_CODE_FREE_MATRIX( descAD );
    PASTE_CODE_FREE_MATRIX( descAUV );
    PASTE_CODE_FREE_MATRIX( descArk );

    return hres;
}


/*------------------------------------------------------------------------
 *  Check the factorization of the matrix A2
 */
static int check_factorization(int N, double *A1, double *A2, int LDA, int uplo, double eps)
{
    double Anorm, Rnorm;
    double alpha;
    int info_factorization;
    int i,j;

    double *Residual = (double *)malloc(N*N*sizeof(double));
    double *L1       = (double *)malloc(N*N*sizeof(double));
    double *L2       = (double *)malloc(N*N*sizeof(double));
    double *work              = (double *)malloc(N*sizeof(double));

    memset((void*)L1, 0, N*N*sizeof(double));
    memset((void*)L2, 0, N*N*sizeof(double));

    alpha= 1.0;

    LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,' ', N, N, A1, LDA, Residual, N);

    /* Dealing with L'L or U'U  */
    if (uplo == MorseUpper){
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'u', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasLeft, CblasUpper, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }
    else{
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L1, N);
        LAPACKE_dlacpy_work(LAPACK_COL_MAJOR,'l', N, N, A2, LDA, L2, N);
        cblas_dtrmm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, N, N, (alpha), L1, N, L2, N);
    }

    /* Compute the Residual || A -L'L|| */
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
           Residual[j*N+i] = L2[j*N+i] - Residual[j*N+i];

    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, N, Residual, N, &Rnorm );
    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, N, A1, LDA, &Anorm );

    printf("============\n");
    printf("Checking the Cholesky Factorization \n");
    printf("-- ||L'L-A||_oo/(||A||_oo.N.eps) = %e eps:%.1e\n",Rnorm/(Anorm*N*eps), eps);

    if ( isnan(Rnorm/(Anorm*N*eps)) || isinf(Rnorm/(Anorm*N*eps)) || (Rnorm/(Anorm*N*eps) > 60.0) ){
        printf("-- Factorization is suspicious ! \n");
        info_factorization = 1;
    }
    else{
        printf("-- Factorization is CORRECT ! \n");
        info_factorization = 0;
    }

    free(Residual); free(L1); free(L2); free(work);

    return info_factorization;
}


/*------------------------------------------------------------------------
 *  Check the accuracy of the solution of the linear system
 */

static int check_solution(int N, int NRHS, double *A1, int LDA, double *B1, double *B2, int LDB, double eps )
{
    int info_solution;
    double Rnorm, Anorm, Xnorm, Bnorm, result;
    double alpha, beta;
    double *work = (double *)malloc(N*sizeof(double));

    alpha = 1.0;
    beta  = -1.0;

    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, NRHS, B2, LDB, &Xnorm );
    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, N,    A1, LDA, &Anorm );
    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, NRHS, B1, LDB, &Bnorm );

    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, N, NRHS, N, (alpha), A1, LDA, B2, LDB, (beta), B1, LDB);
    BLAS_dge_norm( blas_colmajor, blas_inf_norm, N, NRHS, B1, LDB, &Rnorm );

    if (getenv("MORSE_TESTING_VERBOSE")) {
      printf( "||A||_oo=%f\n||X||_oo=%f\n||B||_oo=%f\n||A X - B||_oo=%e eps:%.2e\n", Anorm, Xnorm, Bnorm, Rnorm, eps );
    }
    result = Rnorm / ( (Anorm*Xnorm+Bnorm)*N*eps ) ;
    printf("============\n");
    printf("Checking the Residual of the solution \n");
    printf("-- ||Ax-B||_oo/((||A||_oo||x||_oo+||B||_oo).N.eps) = %e \n", result);

    if (  isnan(Xnorm) || isinf(Xnorm) || isnan(result) || isinf(result) || (result > 60.0) ) {
        printf("-- The solution is suspicious ! \n");
        info_solution = 1;
     }
    else{
        printf("-- The solution is CORRECT ! \n");
        info_solution = 0;
    }

    free(work);

    return info_solution;
}
