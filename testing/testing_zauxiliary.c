/**
 *
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2014 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 *
 **/

/**
 *
 * @file testing_zauxiliary.c
 *
 *  MORSE testing routines
 *  MORSE is a software package provided by Univ. of Tennessee,
 *  Univ. of California Berkeley and Univ. of Colorado Denver
 *
 * @version 0.9.0
 * @author Mathieu Faverge
 * @author Cédric Castagnède
 * @date 2018-11-08
 * @precisions normal z -> c d s
 *
 **/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#if defined( _WIN32 ) || defined( _WIN64 )
#include <windows.h>
#else  /* Non-Windows */
#include <unistd.h>
#include <sys/resource.h>
#endif
#include <morse.h>
#undef REAL
#define COMPLEX
#undef SINGLE
#define DOUBLE
#include "testing_zauxiliary.h"

int   IONE     = 1;
int   ISEED[4] = {0,0,0,1};   /* initial seed for zlarnv() */

int   format[6]= { MorseCM, MorseRM, MorseCCRB, MorseCRRB, MorseRCRB, MorseRRRB };
int   side[2]  = { MorseLeft,    MorseRight };
int   uplo[2]  = { MorseUpper,   MorseLower };
int   diag[2]  = { MorseNonUnit, MorseUnit  };
int   trans[3] = { MorseNoTrans, MorseTrans, MorseConjTrans };
int   itype[3]  = { 1, 2, 3 };
int   storev[2] = { MorseRowwise, MorseColumnwise };
int   norm[4]   = { MorseMaxNorm, MorseOneNorm, MorseInfNorm, MorseFrobeniusNorm };

char *formatstr[6]= { "CM", "RM", "CCRB", "CRRB", "RCRB", "RRRB"};
char *sidestr[2]  = { "Left ", "Right" };
char *uplostr[2]  = { "Upper", "Lower" };
char *diagstr[2]  = { "NonUnit", "Unit   " };
char *transstr[3] = { "N", "T", "H" };
char *itypestr[3]  = { "inv(U')xAxinv(U) or inv(L)xAxinv(L')", "UxAxU' or L'xAxL", "UxAxU' or L'xAxL" };
char *storevstr[2] = { "Rowwise", "Columnwise" };
char *normstr[4]  = { "Max", "One", "Inf", "Fro" };

#define map_cm(m, n, i, j)   ((i) + (j) * (m))
#define map_rm(m, n, i, j)   ((i) * (n) + (j))

int map_CM(int m, int n, int mb, int nb, int i, int j)
{
    int hres = map_cm(m, n, i, j);
    (void)mb;
    (void)nb;
    (void)n;
    return hres;
}

int map_RM(int m, int n, int mb, int nb, int i, int j)
{
    int hres = map_rm(m, n, i, j);
    (void)mb;
    (void)nb;
    (void)m;
    return hres;
}

int map_CCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_CRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_cm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n0 + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m*n0  + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RCRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_cm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_cm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_cm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_cm( m%mb, n%nb, i%mb, j%nb) );
}

int map_RRRB(int m, int n, int mb, int nb, int i, int j) {
    int m0 = m - m%mb;
    int n0 = n - n%nb;
    if ( j < n0 )
        if (i < m0)
            /* Case in A11 */
            return ( map_rm( m/mb, n/nb, i/mb, j/nb )*mb*nb + map_rm( mb,   nb,   i%mb, j%nb) );
        else
            /* Case in A21 */
            return ( m0*n  + ( (j/nb) * (nb*(m%mb)) )       + map_rm( m%mb, nb,   i%mb, j%nb) );
    else
        if (i < m0)
            /* Case in A12 */
            return ( m0*n0 + ( (i/mb) * (mb*(n%nb)) )       + map_rm( mb,   n%nb, i%mb, j%nb) );
        else
            /* Case in A22 */
            return ( m*n0  + (n-n0)*m0                      + map_rm( m%mb, n%nb, i%mb, j%nb) );
}

int (*formatmap[6])(int, int, int, int, int, int) = {  map_CM, map_RM, map_CCRB, map_CRRB, map_RCRB, map_RRRB };

int main (int argc, char **argv)
{
    int ncores, ngpus;
    int info = 0;
    char func[32];

    /* Check for number of arguments*/
    if ( argc < 4) {
        printf(" Proper Usage is : ./ztesting ncores ngpus FUNC ...\n"
               "   - ncores : number of cores\n"
               "   - ngpus  : number of GPUs\n"
               "   - FUNC   : name of function to test\n"
               "   - ... plus arguments depending on the testing function \n");
        exit(1);
    }

    sscanf( argv[1], "%d",   &ncores );
    sscanf( argv[2], "%d",   &ngpus  );
    sscanf( argv[3], "%31s",  func   );

    /* Initialize MORSE */
    /*if(nthreads_per_worker)
     MORSE_InitPar(ncores/nthreads_per_worker, ncudas, nthreads_per_worker);
     else*/
    MORSE_Init( ncores, ngpus);
    MORSE_Disable(MORSE_AUTOTUNING);
    MORSE_Set(MORSE_TILE_SIZE,         32 );
    MORSE_Set(MORSE_INNER_BLOCK_SIZE,   5 );

    argc -= 4;
    argv += 4;
    info  = 0;

    /*
     * Linear system
     */
    if ( strcmp(func, "posv") == 0 ) {
        info += testing_dposv( argc, argv );
    }
    else if ( strcmp(func, "trsmd") == 0 ) {
        info += testing_dtrsmd( argc, argv );
    }
    else {
        fprintf(stderr, "Function unknown\n");
    }

    if ( info == -1 ) {
        printf( "TESTING %s FAILED : incorrect number of arguments\n", func);
    } else if ( info == -2 ) {
        printf( "TESTING %s FAILED : not enough memory\n", func);
    }

    MORSE_Finalize();

    return info;
}
