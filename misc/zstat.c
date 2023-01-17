/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file zstat.c
 *
 * This file contains the function for getting statistics (avg, min, max) of numbers in a matrix.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/

#include <stdio.h>
#include <hicma_common.h>

void zget_stat(HICMA_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat)
{
    double final_avgrank;
    int final_maxrank = 0;
    int minrank = 10000;
    int final_totalrank = 0;
    double *MAT = Ark;
    int64_t i, j, imt, jnt, nelm = 0;
    int ntiles = 0;
    for(imt=0;imt<m;imt++){
        for(jnt=0;jnt<n;jnt++){
            if((uplo == HicmaLower || uplo == HicmaUpper) && imt == jnt)
                continue;
            if((uplo == HicmaUpperLower) && imt == jnt)
                continue;
            if(uplo == HicmaLower && imt < jnt)
                continue;
            if(uplo == HicmaUpper && imt > jnt)
                continue;
            double *A = MAT+imt+jnt*ld;
            int rank = A[0];
            if(rank > final_maxrank){
                final_maxrank = rank;
            }
            if(rank < minrank){
                minrank = rank;
            }
            final_totalrank += rank;
            ntiles++;
            
            if(0){
                //if(jnt<imt)
                printf("%s %d: Tile %d,%d %d ntiles:%d final_totalrank:%d\n", __FILE__, __LINE__, imt, jnt, rank, ntiles, final_totalrank);
            }
        }
    }
    final_avgrank = (final_totalrank)/(ntiles*1.0);
    stat->min = minrank;
    stat->max = final_maxrank;
    stat->avg = final_avgrank;
}
void zprint_stat(HICMA_stat_t stat)
{
    printf("avg:%g min:%d max:%d\n", stat.avg, stat.min, stat.max);
}

