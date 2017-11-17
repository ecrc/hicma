/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file descprint.c
 *
 * This file contains the functions for printing numerical values stored inside MORSE descriptors. 
 * 
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/
#include "morse.h"
#include "auxdescutil.h"
int64_t nelm_limit = 100000;
int64_t nrow_limit = 100;
int64_t ncol_limit = 100;
void printmat(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics){
    printf("M:%d N:%d LD:%d\n", m, n, ld);
    int64_t i, j, nelm = 0;
    for(i=0;i<m;i++){
        int ncols = 0;
        for(j=0;j<n;j++){
            //if(A[j*ld+i] > 1e-6)
            //printf(" (%d,%d):%+.2e", i,j,A[j*ld+i]);
            printf("%+.2e", A[j*ld+i]);
            //printf("%g ", A[j*tld(descZ)+i]);
            //printf("%g\t", A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
            nelm++;
            if(nelm >= nelm_limit){
                printf("\n");
                return;
            }
            if(j==ncol_limit)
                break;
            if((j+1)%ics==0)
                printf("|");
            ncols+= 9;
        }
        printf("\n");
        if(i==nrow_limit)
            break;
        if((i+1)%irs==0){
          for(j = 0; j < ncols; j++){
              printf("-");
          }
          printf("\n");
        }
    }
}
void printdescrk(MORSE_desc_t *descZ, int64_t rank){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            printf("Tile %d,%d\n", imt, jnt);
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printmat(A, descZ->nb, rank, tld(descZ), -1, -1);
        }
    }
}
void printdesc(MORSE_desc_t *descZ){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            printf("Tile %d,%d\n", imt, jnt);
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printmat(A, descZ->nb, descZ->nb, tld(descZ), -1, -1);
        }
    }
}
