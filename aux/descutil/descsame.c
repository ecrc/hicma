/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file descsame.c
 *
 * This file contains the functions for checking whether two MORSE descriptors are same or not.
 * 
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#include "morse.h"
#include "auxdescutil.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
/** 
 * This function isused to compare to descriptors in terms of size and numerical values.
 * Do not use this function with MPI
 */
void check_same(MORSE_desc_t *descL, MORSE_desc_t *descR, char diag, char uplo){
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
            if(uplo == 'L' && imt < jnt){
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
void check_same_array(double *L, double *R, int nelm, int line, char *file){
    int i;
    int error_encountered = 0;
    for(i=0; i<nelm; i++){
        double valL = L[i];
        double valR = R[i];
        double diff = fabs(valL - valR);
        double thresh = 1e-12;
        //double thresh = 1e-2;
        if(diff > thresh ){
            printf("Elm:%d diff:%g val:%.14e %.14e\t",  i,  diff, valL, valR);
            error_encountered = 1;
            //exit(1);
            break;
        }
    }
    if (error_encountered == 0) {
        printf("arrays are same at line %d of %s\n", line, file);
    } else {
        printf("arrays are NOT SAME !!!!!!!!!!!!!!!!!!!!!!  at line %d of %s\n", line, file);
        exit(1);
    }
}
