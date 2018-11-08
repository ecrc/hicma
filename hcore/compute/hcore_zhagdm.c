/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hcore_zhagdm.c
 *
 *  HiCMA HCORE routine for generating dense matrix.
 *  
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
//#include "hcore/include/hcore.h"
#include "morse.h"
#include "hcore_z.h"
#include <assert.h>
#include <stdio.h>
#include <sys/time.h>//FIXME for gettimeofday

#include "hicma.h"
#include "hicma_common.h"

#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-randtlr.h"
#ifdef LAPACKE_UTILS
#include <lapacke_utils.h>
#endif
#include "coreblas/coreblas.h"
#include "coreblas/lapacke.h"

extern void _printmat(double * A, int m, int n, int ld);

void HCORE_zhagdm( 
        int nrows_Dense,
        int ncols_Dense,
        double *Dense,
        int ld_Dense,
        int tile_row_index,
        int tile_col_index
        )
{
    struct timeval tvalBefore, tvalAfter; 
    gettimeofday (&tvalBefore, NULL);
    if(HICMA_get_print_index() == 1){
        fprintf(stderr, "%d+HAGDM\t|(%d,%d) m:%d n:%d ld:%d\n",MORSE_My_Mpi_Rank(), 
                tile_row_index, tile_col_index, 
                nrows_Dense, ncols_Dense, ld_Dense);
    }
    STARSH_cluster *RC = HICMA_get_starsh_format()->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;
    HICMA_get_starsh_format()->problem->kernel(nrows_Dense, ncols_Dense, 
            RC->pivot+RC->start[tile_row_index], 
            CC->pivot+CC->start[tile_col_index],
            RD, CD, Dense, ld_Dense);
    if(HICMA_get_print_mat() == 1){
        printf("%d\tHAGDM-DENSE-output\n", __LINE__);
        _printmat(Dense, nrows_Dense, ncols_Dense, ld_Dense);
    }
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-HAGDM\t|(%d,%d) m:%d n:%d ld:%d\t\t\t\t\tHAGDM: %.4f\n",
                MORSE_My_Mpi_Rank(), 
                tile_row_index, tile_col_index, 
                nrows_Dense, ncols_Dense, ld_Dense,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
               );
    }
}



