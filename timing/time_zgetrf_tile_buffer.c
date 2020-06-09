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

int print_progress = 0;   // Print progress about the execution
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
void AL4SAN_near_sca(double _Complex *rhsresults, int *nip, int *ntrian);

#if defined(HICMA_COMPLEX)
    int
RunTest(int *iparam, double *dparam, morse_time_t *t_, char* rankfile)
{
  PASTE_CODE_IPARAM_LOCALS( iparam );  
  //HICMA_set_print_mat();
    /*HICMA_set_print_index();*/
   printf("\n ************************I am in complex part ************************\n");
   PROGRESS("RunTest started");

   // int nip = 3;
   // int ntrian = 120;

    //int nip = NB;
    //int ntrian = N/NB;
    //int ntrian = 6; //obtained from command line arguments

    printf("\n nip:%d , ntrian:%d N:%d, NB:%d\n", nip, ntrian, N, NB);
//ka_ex
  //allocating complex matrix

   //printf ( "Values are: %g%+gi %g%+gi  %d\n", creal(zzsparse_tempt[0]), cimag(zzsparse_tempt[0]),creal(zzsparse_tempt[1]), cimag(zzsparse_tempt[1]));

    //AL4SAN_Matrix_Generation(&nip, &ntrian, zzsparse_tempt, crhs);
    MORSE_user_tag_size(31,27);

    //chameleon/runtime/starpu/control/runtime_descriptor.c
    //MORSE_user_tag_size(31,26);
    //MORSE_user_tag_size(31,29);
    int saveNB_2 = NB;
    NB = MB;

   int global_nt=N/NB;
    int local_nt=NB/nip;
   printf("\n global_nt:%d local_nt:%d N:%d NB:%d \n", global_nt, local_nt, N, NB);
   Acoustic_Init(&nip, &ntrian);
   //printf ( "Hala Values are: %g%+gi %g%+gi \n", creal(zzsparse_tempt[0]), cimag(zzsparse_tempt[0]),creal(zzsparse_tempt[1]), cimag(zzsparse_tempt[1]));

   PASTE_CODE_ALLOCATE_MATRIX_TILE( descDense, solve, MORSE_Complex64_t, MorseComplexDouble, (nip* ntrian), (nip* ntrian), (nip* ntrian) ); 
    PROGRESS("descDense is allocated");

    PROGRESS("HICMA_zgene started");
    //HICMA_zgene_Tile( ntrian,  nip,  local_nt,  NB, descDense);
    if(0)printf("%s %d: descACO LDA:%d M:%d MB:%d NB:%d\n", __FILE__, __LINE__, (nip* ntrian), (nip* ntrian), MB, NB); fflush(stdout);
    NB = saveNB_2;
    //printf("\n descACO->nt:%d, NB:%d\n", descACO->nt, NB);
    // print progress info only on ROOT process
    if(MORSE_My_Mpi_Rank() != 0)
        print_progress = 0;
    PROGRESS("HICMA_zgene ended");

    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);

    // this paramater enables storing only diagonal tiles in a tall and skinny matrix
    store_only_diagonal_tiles = 1;
    //chameleon/runtime/starpu/control/runtime_descriptor.c
    //MORSE_user_tag_size(31,26);
    //MORSE_user_tag_size(31,29);
//    MORSE_user_tag_size(31,27);// When I added tile_to_lapack for descArk, I got not enough number of desc error

    //MORSE_user_tag_size(31,27);// When I added tile_to_lapack for descArk, I got not enough number of desc error
    // get parameters coming from command line
//    PASTE_CODE_IPARAM_LOCALS( iparam );

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
    //printf("\n N = ncols_AD:%d\n", N);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAD, 1, MORSE_Complex64_t, MorseComplexDouble, LDA, M, N );
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
      
    //complex desc
    //PASTE_CODE_ALLOCATE_MATRIX_TILE( descDense, 1, MORSE_Complex64_t, MorseComplexDouble, (nip* ntrian), (nip* ntrian), (nip* ntrian) );
    //MORSE_Lapack_to_Tile(zzsparse_tempt, (nip* ntrian), descACO);



     //uncommented if you want to print the tiles
   //PASTE_TILE_TO_LAPACK(descACO, ACO, 1, MORSE_Complex64_t, (nip* ntrian), (nip* ntrian));
    //_printmat_complex(ACO, (nip* ntrian), (nip* ntrian), (nip* ntrian));
    //free(ACO);

    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);
    //PASTE_CODE_ALLOCATE_MATRIX_TILE( descDense, 1, MORSE_Complex64_t, MorseComplexDouble, ld_Dense, ncols_Dense, ncols_Dense );

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
    //printf("N:%d NB:%d\n", N, NB);
    saveN = N;
    N = N * 2;
    saveNB = NB;
    NB = NB * 2;
    //printf("N:%d NB:%d\n", N, NB);
    //PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUV, 1, MORSE_Complex64_t, MorseComplexDouble, ld_AUV, nrows_AUV, N );
    N = saveN;
    NB = saveNB;

    /* tile dimension of rank descriptor must be 1 */
    /* when LD for rk matrices is 1, program exits*/
    int bigMB = MB;
    int bigNB = NB;
    MB = NB = 1;
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descArk, 1, double, MorseRealDouble, MT, MT, NT);
    PROGRESS("descArk is allocated");
    MB = bigMB;
    NB = bigNB;

    int diag_dense = 1;
    int fixedrank = iparam[IPARAM_RK]; //genargs->k
    double fixedacc = pow(10, -1.0*iparam[IPARAM_ACC]);

    char sym;
    if (run_getrf)
        sym = 'S';
    else
        sym = 'N';
    int probtype = iparam[IPARAM_HICMA_STARSH_PROB];
    int maxrank  = iparam[IPARAM_HICMA_STARSH_MAXRANK];
    //double ddecay = pow(10, -1.0*iparam[IPARAM_HICMA_STARSH_DECAY]);
    double ddecay = dparam[IPARAM_HICMA_STARSH_DECAY];


    int  initial_maxrank, final_maxrank;
    double initial_avgrank, final_avgrank;
    struct timeval tvalBefore, tvalAfter;  // removed comma

    if(0)printf("check N:%d, NB:%d", N, NB); 
    //N=M;
    //NB=MB;
   int initmaxrank=NB;
    if(0)printf("\n %s %d maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d\n", __FILE__, __LINE__,  initmaxrank, N, NB, M, MB, ld_AUV, NT, MT);
    int isPerused=iparam[IPARAM_HICMA_ISPRECENT];
    precent1=((float)iparam[IPARAM_HICMA_PRECENT1])/100;
    precent2=((float)iparam[IPARAM_HICMA_PRECENT2])/100;
    int sp =precent1*initmaxrank; if (sp==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); sp=1;}
    //int ts=2*(initmaxrank+sp);
    int ts=(initmaxrank+sp);
    if(0)printf("\n computed ts:%d", ts);
    ts=chameleon_min(ts, 2*NB);
    //int as=(2*(initmaxrank+sp))*NT;
    int as=(initmaxrank+sp)*NT;
    if(0)printf(" and computed as:%d\n", as);
    as=chameleon_min(2*N, as);
    if(0) printf("\n ts:%d, NB:%d, as:%d, N:%d\n", ts, NB, as, N);
    if(chameleon_min(ts, 2*NB)==2*NB || chameleon_min(2*N, as)==2*N || isPerused==0)
      {
          printf("percentage of inital max rank is not used since it is bigger than 2*MB\n"); 
          isPerused=0;
          NB=2*NB;
          N=2*N;
     } 
     else
      {
        NB=ts;
	N=as;
      }
    isPrecentageused=isPerused;

    if(0)printf("\n %s %d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d isPerused:%d\n", __FILE__, __LINE__, N, NB, M, MB, ld_AUV, NT, MT, isPerused);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUV, 1, MORSE_Complex64_t, MorseComplexDouble, ld_AUV, nrows_AUV, N );
    PROGRESS("descAUV is allocated");
    if(0)printf("\n %s %d  maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d descAUV->lmt:%d, descAUV->lnt:%d descAUV->mt:%d, descAUV->nt:%d\n", __FILE__, __LINE__,  maxrank, N, NB, M, MB, ld_AUV, NT, MT, descAUV->lmt,  descAUV->lnt, descAUV->mt,  descAUV->nt);
   descAUV->lmt=MT;
   descAUV->lnt=NT;
   descAUV->mt=MT;
   descAUV->nt=NT;
    if(0)printf("\n %s %d  maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d descAUV->lmt:%d, descAUV->lnt:%d descAUV->mt:%d, descAUV->nt:%d\n", __FILE__, __LINE__,  maxrank, N, NB, M, MB, ld_AUV, NT, MT, descAUV->lmt,  descAUV->lnt, descAUV->mt,  descAUV->nt);
    N = saveN;
    NB = saveNB;
    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);

	// DO NOT enforce compression of diagonal tiles
    int compress_diag = 0;
    PROGRESS("nompi zgytlr starting");
    //descDense original problem
    gettimeofday (&tvalBefore, NULL);
    HICMA_zgytlr_Tile_buffer(MorseUpperLower, descAUV, descAD, descArk, 0, maxrank, fixedacc, compress_diag, ntrian,  nip,  local_nt,  MB);//, descDense);
    gettimeofday (&tvalAfter, NULL);
   // printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);
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
       int newMaxRank;
    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descArk, Ark_initial, 1, double, MT, NT );
        if(MORSE_My_Mpi_Rank()==0){

            sprintf(rankfile, "%s_initialranks", rankfile);
            fwrite_array(descArk->m, descArk->n, descArk->m, Ark_initial, rankfile);
            print_array(descArk->m, descArk->n, descArk->m, Ark_initial, stdout);

            HICMA_stat_t hicma_statrk_initial;
            zget_stat(MorseUpperLower, Ark_initial, MT, NT, MT,  &hicma_statrk_initial);
            if(1)printf("initial_ranks:");
            if(1)zprint_stat(hicma_statrk_initial);
             newMaxRank=hicma_statrk_initial.max;
            fflush(stderr);
            fflush(stdout);
        }
   }

    //N=M;
    //NB=MB;
    //printf("\n %s %d 2maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d\n", __FILE__, __LINE__,  newMaxRank, N, NB, M, MB, ld_AUV, NT, MT);
    //isPerused=iparam[IPARAM_HICMA_ISPRECENT];
    //precent2=((float)iparam[IPARAM_HICMA_PRECENT2])/100;
    float drivedper=(float)newMaxRank/(float)MB;
     //printf("\n drivedper:%f\n", drivedper);
    sp =newMaxRank*drivedper; if (sp==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); sp=1;}
    ts=2*(newMaxRank+sp);
    ts=chameleon_min(ts, 2*NB);
    as=(2*(newMaxRank+sp))*NT;
    as=chameleon_min(2*N, as);
     if(0) printf("\n ts:%d, NB:%d, as:%d, N:%d\n", ts, NB, as, N);
    if(chameleon_min(ts, 2*NB)==2*NB || chameleon_min(2*N, as)==2*N || isPerused==0)
      {
          printf("percentage of new max rank is not used since it is bigger than 2*MB\n");
          isPerused=0;
          NB=2*NB;
          N=2*M;
      }
     else
     {
      NB=ts;
      N=as;
     }
    isPrecentageused=isPerused;

if (isPrecentageused){
    if(0)printf("\n %s %d 2N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d isPrecentageused:%d\n", __FILE__, __LINE__, N, NB, M, MB, ld_AUV, NT, MT, isPrecentageused);
    PASTE_CODE_ALLOCATE_MATRIX_TILE( descAUVnew, 1, MORSE_Complex64_t, MorseComplexDouble, ld_AUV, nrows_AUV, N );
    if(0)printf("\n %s %d  2maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d descAUV->lmt:%d, descAUV->lnt:%d descAUV->mt:%d, descAUV->nt:%d\n", __FILE__, __LINE__,  newMaxRank, N, NB, M, MB, ld_AUV, NT, MT, descAUVnew->lmt,  descAUVnew->lnt, descAUVnew->mt,  descAUVnew->nt);
   descAUVnew->lmt=MT;
   descAUVnew->lnt=NT;
   descAUVnew->mt=MT;
   descAUVnew->nt=NT;
    if(0)printf("\n %s %d  2maxrank:%d N:%d, NB:%d M:%d MB:%d ld_AUV:%d NT:%d, MT:%d descAUV->lmt:%d, descAUV->lnt:%d descAUV->mt:%d, descAUV->nt:%d\n", __FILE__, __LINE__,  newMaxRank, N, NB, M, MB, ld_AUV, NT, MT, descAUVnew->lmt,  descAUVnew->lnt, descAUVnew->mt,  descAUVnew->nt);
    N = saveN;
    NB = saveNB;
    PROGRESS("nompi zcpy starting");
    //descDense original problem
    gettimeofday (&tvalBefore, NULL);
    HICMA_zcpy_Tile(MorseUpperLower, descAUV, descAUVnew, descArk, 0, newMaxRank, fixedacc, compress_diag);//, isPerused);
    PASTE_CODE_FREE_MATRIX( descAUV );
    gettimeofday (&tvalAfter, NULL);
    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);
    if(MORSE_My_Mpi_Rank()==0){
        printf("Tcpy:%g\n",
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
        fflush(stderr);
        fflush(stdout);
    }
    PROGRESS(" zcpy finished");
  }
//  else { printf("\n descAUVnew=descAUV\n");descAUVnew=descAUV;}


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
    //descDense = descACO;
    if(0)printf("%s %d: descDense LDA:%d M:%d MB:%d NB:%d\n", __FILE__, __LINE__, LDA, M, MB, NB); fflush(stdout);
    //PASTE_TILE_TO_LAPACK( descDense, Adense, check, MORSE_Complex64_t, LDA, M );
    //if(main_print_mat){printf("Adense\n");_printmat_complex(Adense,M,M,LDA);}
    MORSE_Complex64_t one = 1.0, zero = 0.0, minusone = -1.0, diagVal = M;
    MORSE_Complex64_t* swork = NULL;
    //double* cp_L_Adense =  calloc(LDA*M, sizeof(double)); 
    /*if(check && P*Q==1){
        swork  = calloc(2*M, sizeof(MORSE_Complex64_t));
        {size_t i, j;
            MORSE_Complex64_t* orgAdense = calloc(LDA*M, sizeof(MORSE_Complex64_t));
            for(j = 0; j < M; j++){
                for(i = 0; i < M; i++){
                    orgAdense[j*LDA+i] = Adense[j*LDA+i];
                }
            }
//printmat(orgAdense,M,M,LDA,MB, MB);
int *ipiv=(int*)malloc(chameleon_max(1, chameleon_min(M, M))*sizeof(int));
            int info = LAPACKE_dgetrf_work(
                    LAPACK_COL_MAJOR,M,
                    M, orgAdense, LDA,ipiv);
 free(ipiv);
//printf("\n checkcing after LU gerf \n");
//printmat(orgAdense,M,M,LDA,MB, MB);
//exit(0);
            if(info != 0){ 
                fprintf(stderr, "%s\t|%d\t|Error in LAPACK getrf. info:%d, This errors means "
                        "that the matrix generated is not positive definite\n", __FILE__, __LINE__, info);
            }
            //for(j = 0; j < M; j++){
              //  for(i = 0; i < j; i++){
                //    orgAdense[j*LDA+i] = zero;
                //}
            //}

            //for(j = 0; j < M; j++) {                
                //for(i = 0; i < M; i++){
                    //cp_L_Adense[j*LDA+i] = orgAdense[j*LDA+i];
                //}
            //}
            if(main_print_mat ){printf("L of Adense\n");printmat(orgAdense,M,M,LDA,MB, MB);}
             MORSE_Complex64_t normOrgAdense = 0.0;
            //HICMA_znormest(M, M, orgAdense, &normOrgAdense, swork);
            //printf("norm_L_OrgAdense:%e\n",normOrgAdense);
            free(orgAdense);
        }
    }*/
    //PASTE_TILE_TO_LAPACK( descDense, Adense2, check, MORSE_Complex64_t, LDA, M );
    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);
    PROGRESS("pasting original dense descAD into Adense and Adense2 finished");
    PROGRESS("getrf started");
   // printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    //exit(-1);
    START_TIMING();
    HICMA_zgetrf_Tile(
            MorseLower, /* @noha: this paramater must be removed*/
            descAUV, descAD, descArk, fixedrank, maxrank, fixedacc); //, isPerused);
    STOP_TIMING();
    fflush(stderr);
    fflush(stdout);
    //PROGRESS("getrf finished");
   // printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
   // fflush(stdout);
    // exit(-1);
//printf("%s %d\n", __FILE__, __LINE__);exit(0);
 /*   if(check){
        HICMA_zuncompress(MorseUpperLower, descAUV, descDense, descArk);
        HICMA_zdiag_vec2mat(descAD, descDense);
        PASTE_CODE_FREE_MATRIX( descAD ); //@KADIRLBL001  
        descAD = descDense; // descAD was only diagonals.
        // After this line, descAD is dense matrix containing approximate L
        // So no need to adapt below code for descAD containg only diagonals.
    }*/
    //PROGRESS("getrf finished");
    //printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
    //fflush(stdout);
    if(calc_rank_stat == 1) {
        PASTE_TILE_TO_LAPACK( descArk, Ark_final, 1, double, MT, NT );
        if(MORSE_My_Mpi_Rank()==0){
            sprintf(rankfile, "%s_finalranks", rankfile);
            fwrite_array(descArk->m, descArk->n, descArk->m, Ark_final, rankfile);
            print_array(descArk->m, descArk->n, descArk->m, Ark_final, stdout);
            HICMA_stat_t hicma_statrk_final;
            zget_stat(MorseUpperLower, Ark_final, MT, NT, MT,  &hicma_statrk_final);
            if(1)printf("final_ranks:");
            if(1)zprint_stat(hicma_statrk_final);
            if(hicma_statrk_final.max > MB) {
                printf("Maximum actual rank %d is larger than block size %d so this program won't give correct results. Exiting...\n", hicma_statrk_final.max, MB);
                exit(-1);
            }
            fflush(stderr);
            fflush(stdout);
        }
    }

    int check_dense = 0;
    int check_app = 1;
    if(solve){
        HICMA_zuncompress(MorseUpperLower, descAUV, descDense, descArk);
        HICMA_zdiag_vec2mat(descAD, descDense);
        PASTE_CODE_FREE_MATRIX( descAD ); //@KADIRLBL001
        descAD = descDense; // descAD was only diagonals.
        // After this line, descAD is dense matrix containing approximate L
        // So no need to adapt below code for descAD containg only diagonals.
    }

   if (solve){
       //PASTE_TILE_TO_LAPACK( descDense, checklu, 1, MORSE_Complex64_t, (nip* ntrian), (nip* ntrian) );
       //_printmat_complex(checklu,  (nip* ntrian), (nip* ntrian), (nip* ntrian));
//       double _Complex *crhs;
//       crhs=malloc((nip* ntrian)*sizeof (double _Complex ));
//       AL4SAN_RHS( crhs , &ntrian, &nip);
       //_printmat_complex(crhs,  (nip* ntrian), 1, (nip* ntrian));
      NB=MB;
       PASTE_CODE_ALLOCATE_MATRIX_TILE( descRHS, 1, MORSE_Complex64_t, MorseComplexDouble, (nip* ntrian), (nip* ntrian), 1 );
       //MORSE_Lapack_to_Tile(crhs, (nip* ntrian) , descRHS);
      if(0)printf("\nntrian:%d,  nip:%d,  local_nt:%d, MB:%d\n",ntrian,  nip,  local_nt, MB);
      AL4SAN_RHS_zgene_Tile( ntrian,  nip,  local_nt, MB, descRHS);
       MORSE_ztrsm_Tile(MorseLeft, MorseLower, MorseNoTrans, MorseUnit, (MORSE_Complex64_t)1.0, descDense, descRHS);
       MORSE_ztrsm_Tile(MorseLeft, MorseUpper, MorseNoTrans, MorseNonUnit, (MORSE_Complex64_t)1.0, descDense, descRHS);
       PASTE_TILE_TO_LAPACK( descRHS, resultsrhs, 1, MORSE_Complex64_t, (nip* ntrian), 1);
       if(0)_printmat_complex(resultsrhs,  (nip* ntrian), 1, (nip* ntrian));
       PASTE_CODE_FREE_MATRIX( descRHS );
       FILE *f;
       f = fopen("rhs.bin","wb");
       fwrite(resultsrhs,sizeof(double _Complex),(nip* ntrian),f);
       AL4SAN_near_sca(resultsrhs, &ntrian, &nip);
       free(resultsrhs);
     }

/*    if(check == 0){
        check_dense = check_app = 0;
    }
  if(check_app ) {
        PROGRESS("checking accuracy");
        if( MORSE_My_Mpi_Rank()==0){
#ifndef COMPLEX
            if(main_print_mat){printf("Adense2\n");printmat(Adense2,M,M,LDA,MB, MB);}
            double normA;
            PROGRESS("normaA started");
            //HICMA_znormest(M, M, Adense2, &normA, swork);
            normA = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, M, M, Adense2, M);
            // Ahicma: result of TLR getrf
            PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, MORSE_Complex64_t, LDA, M );
            //if(0){size_t i,j;
                //for(j = 0; j < M; j++) {               
                    //for(i = 0; i < M; i++){
                        //Ahicma[j*LDA+i] = cp_L_Adense[j*LDA+i];
                    //}
                //}
            //}
            double normAhicma = 0.0;
            {size_t i, j;
                MORSE_Complex64_t* orgAhicma = calloc(LDA*M, sizeof(MORSE_Complex64_t));
                for(j = 0; j < M; j++){
                    for(i = 0; i < M; i++){
                        orgAhicma[j*LDA+i] = Ahicma[j*LDA+i];
                    }
                }
                //HICMA_znormest(M, M, orgAhicma, &normAhicma, swork);
                normAhicma = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, M, M, orgAhicma, M);
                free(orgAhicma);
            }
            if(set_diag){size_t j; for(j = 0; j < M; j++){ Ahicma[j*LDA+j] = diagVal; } }
//printf("Setting upper part of AhicmaL to zero\n"); fflush(stdout);
            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = 0; i <= j; i++){
                        Ahicma[j*LDA+i] = zero;
                    }
                    Ahicma[j*LDA+j] = one;
                }
            }
            if(main_print_mat){printf("AhicmaL\n");_printmat_complex(Ahicma ,M,M,LDA);} //ow rank L
            //LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', M, Ahicma, LDA);
            // AhicmaT: transpose of Ahicma
            PROGRESS("copy descAd into AhicmaT started");
            PASTE_TILE_TO_LAPACK( descAD,  AhicmaT , check, MORSE_Complex64_t, LDA, M ); //low rank U


            {size_t i, j;
                for(j = 0; j < M; j++){
                    for(i = j+1; i < M; i++){
                        AhicmaT[j*LDA+i] = zero;
                    }
                }
            }
            //if(main_print_mat){printf("Ahicma-upperzero\n");printmat(Ahicma,M,M,LDA, MB, MB);}
            //LAPACKE_dge_trans(LAPACK_COL_MAJOR, M, M, Ahicma, LDA, AhicmaT, LDA);
            if(main_print_mat){printf("AhicmaU\n");_printmat_complex(AhicmaT,M,M,LDA);}

            PROGRESS("TRMM started");

            cblas_ztrmm (CblasColMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, M, M, &one, Ahicma, LDA, AhicmaT, LDA);
            if(main_print_mat){printf("Ahicma*AhicmaT\n");_printmat_complex(AhicmaT,M,M,LDA);}
            //double tmpnorm;normest(M, M, AhicmaT, &tmpnorm, swork);printf("tmpnorm:%e\n",tmpnorm);

            size_t nelm = M * M;
            if(main_print_mat)printf("nelm:%zu M:%d N:%d\n", nelm, M, N);
            PROGRESS("DAXPY started");
            if(main_print_mat){printf("Adense\n");_printmat_complex(Adense,M,M,LDA);}
            cblas_zaxpy(nelm, &minusone, AhicmaT, 1, Adense, 1);
            if(main_print_mat){printf("Adense-(Ahicma*AhicmaT)\n");_printmat_complex(Adense,M,M,LDA);}

            double normDenseAppDiff;
            PROGRESS("Norm of difference started");
            //HICMA_znormest(M, M, Adense, &normDenseAppDiff, swork);
            //double LAPACKE_zlange (int matrix_layout, char norm, lapack_int m, lapack_int n, const lapack_complex_double * a, lapack_int lda);
            normDenseAppDiff = LAPACKE_zlange(LAPACK_COL_MAJOR, norm, M, M, Adense, M);
            double accuracyDenseAppDiff = normDenseAppDiff/normA;
            printf("normA:%.2e normDenseAppdiff:%.2e Accuracy: %.2e\n", normA, normDenseAppDiff,  accuracyDenseAppDiff);
            dparam[IPARAM_RES] = normDenseAppDiff;
            dparam[IPARAM_ANORM] = normA;
            dparam[IPARAM_XNORM] = normA;
            dparam[IPARAM_BNORM] = normAhicma;
#endif
        } else {
            PASTE_TILE_TO_LAPACK( descAD, Ahicma, check, MORSE_Complex64_t, LDA, M );
            PASTE_TILE_TO_LAPACK( descAD,  AhicmaT, check, MORSE_Complex64_t, LDA, M );
        }
        PROGRESS("checking accuracy is finished");
    }*/

    PASTE_CODE_FREE_MATRIX( descAUV );
    //PASTE_CODE_FREE_MATRIX( descACO); /** Do not free descACO, it has an alias descDense */
    PROGRESS("descAUV is freed");
    if(check == 0 && solve ==0) { // If there is no check, then descAD and descDense are different. Refer to @KADIRLBL001
        PASTE_CODE_FREE_MATRIX( descAD );
        PROGRESS("descAD is freed");
    }
    
    PASTE_CODE_FREE_MATRIX( descArk );
    PROGRESS("descArk is freed");
if(solve) {PASTE_CODE_FREE_MATRIX( descDense );
    PROGRESS("descDense is freed");}
    PROGRESS("freed descs");


    return 0;
}

#endif
