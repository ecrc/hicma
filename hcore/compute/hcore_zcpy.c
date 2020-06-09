/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zgytlr.c
 *
 *  HiCMA HCORE routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 * @precisions normal z -> c d s
 **/
//#include "hcore/include/hcore.h"
#include "morse.h"
#include "hcore_z.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday
#include <stdlib.h>//FIXME for malloc

 #define Stringize( L )     #L 
 #define MakeString( M, L ) M(L)
 #define $Line MakeString( Stringize, __LINE__ )
 #define Reminder __FILE__ "(" $Line ") : Reminder: "

#define COMPLEX
#undef REAL
//  HCORE_zgytlr - Generate a tile for random matrix.

#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"
#include "hicma_common.h"
//extern int isPrecentageused;
//extern float precent2;

extern STARSH_blrf *mpiF;
extern int print_index;
int print_index_end;
extern int store_only_diagonal_tiles;
extern int global_check;
extern int print_mat;



extern void _printmat_complex(MORSE_Complex64_t * A, int m, int n, int ld);

//extern int global_always_fixed_rank; //FIXME this does not work, I do not know why
//extern int global_fixed_rank;
int global_always_fixed_rank;
int global_fixed_rank;

void HCORE_zcpy( int m, int n, /*dimension of squareAD*/
        int nAUV, int nAUVnew, 
        MORSE_Complex64_t *AUV,
        MORSE_Complex64_t *AUVnew,
        int Ark, 
        int ldamUV,
        int ldamUVnew,
        int bigM, int ii, int jj, unsigned long long int seed,
        int maxrank, double tol, int compress_diag
        )
{
    int64_t i, j;
    //printf("m:%d n:%d bigM:%d m0:%d n0:%d\n", m, n, bigM, m0, n0);
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    if(print_index){
        fprintf(stderr, "%d+zcopy\t|(%d,%d) m:%d n:%d ldamUV:%d ldamUV:%d\n",MORSE_My_Mpi_Rank(), ii, jj, m, n,  ldamUV, ldamUVnew);
    }

    int rank = Ark;

    

    if(print_mat){
        printf("%d\told-AUV-input\n", __LINE__);
        //_printmat_complex(AD, m, n, lda);
        //_printmat_complex(AU, m, rank, ldu);
        //_printmat_complex(AV, ldv, rank, ldv);
        _printmat_complex(AUV, m, nAUV, ldamUV);
    }
    char chall = 'A';
    if(ii == jj){ /** copy from Dense to AD */
        return;
    }
    int sp=precent2*rank;
    if (sp==0) {printf("%s %s %d: required percenatge from rank is zero, we will but at least one col\n", __FILE__, __func__, __LINE__); sp=1;}

    int ncols=2*(rank+sp);
    if(print_index){
        printf("\n%s %d rank:%d, sp:%d, precent2:%f, ncols:%d\n", __FILE__, __LINE__, rank, sp, precent2, ncols);
    }

    LAPACK_zlacpy(&chall,
                &m, &ncols, AUV, &ldamUV, AUVnew, &ldamUVnew);


    if(print_mat){
        printf("%d\tnew-AUV-output with rank:%d\n", __LINE__, rank);
        /*_printmat_complex(AD, m, n, lda);*/
        _printmat_complex(AUVnew, m, nAUVnew, ldamUVnew);
    }
    if(print_index || print_index_end ){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-zcpy\t|(%d,%d) rk:%d m:%d n:%d ldamUV:%d ldamUV:%d\t\t\t\t\tzcpy: %.4f\n",MORSE_My_Mpi_Rank(),ii,jj,Ark,m, n, ldamUV, ldamUVnew,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                 +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
                );
    }
    return;
}




