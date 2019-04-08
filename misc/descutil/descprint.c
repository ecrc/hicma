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
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#include "morse.h"
#include "auxdescutil.h"
int64_t nelm_limit = 100000;
int64_t nrow_limit = 100;
int64_t ncol_limit = 100;
void printmat_format(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics, int format){
    printf("%s@%d: M:%d N:%d LD:%d", __FILE__, __LINE__, m, n, ld);
    printf("\n");
    int64_t i, j, nelm = 0;
    for(i=0;i<m;i++){
        if(format == 1 && i==0){
            printf("[");
        }
        if(format == 1){
            printf("[");
        }
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
            if(format == 1 && j != n-1){
                printf(",");
            }
            nelm++;
            if(nelm >= nelm_limit){
                printf("\n");
                return;
            }
            if(j==ncol_limit)
                break;
            if(format == 1){
            } else {
                if((j+1)%ics==0)
                    printf("|");
            }
            ncols+= 9;
        }
        if(format == 1){
            printf("]");
        }
        if(format == 1 && i != m-1 ){
            printf(",");
        }
        if(format == 1 && i == m-1){
            printf("]");
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
void printmat(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics){
    printmat_format(A, m, n, ld, irs, ics, 0);
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
#define tld(d) (d->mb)
#define tsa(d,i,j) (((j)*(d->mt)+(i))*(d->mb)*(d->nb))
#define A(m,n) A,  m,  n
int64_t _nelm_limit = 1200;
int64_t _nrows_limit = 37;
int64_t _ncols_limit = 33;
int64_t _ndim_limit = 37;
#include <stdio.h>
void _printmat(double * A, int64_t m, int64_t n, int64_t ld){
    printf("%s@%d M:%d N:%d LD:%d %p [\n", __FILE__, __LINE__, m, n, ld, A);
    int64_t i, j, nelm = 0;
    for(i=0;i<m;i++){
        printf("[");
        for(j=0;j<n;j++){
            printf("%+.4e", A[j*ld+i]);
            //printf("%g ", A[j*tld(descZ)+i]);
            //printf("%g\t", A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
            //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
            //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
            if(j!=n-1){
                printf(",");
            }
            nelm++;
            if(nelm >= _nelm_limit){
                printf("\n");
                return;
            }
            if(j==_ncols_limit)
                break;
        }
        printf("]");
        if(i!=m-1){
            printf(",");
            printf("\n");
        }
        //printf("\n");
        if(i==_nrows_limit)
            break;
    }
    printf("]\n");
}
void _printdescs(MORSE_desc_t *descD,MORSE_desc_t *descU, MORSE_desc_t *descV,  MORSE_desc_t *descRk){
    int64_t i, j, imt, jnt;
    printf("\n");
    for(imt=0;imt<descD->mt;imt++){
        for(jnt=0;jnt<descD->nt;jnt++){
            if(imt < jnt) continue;
            double *MAT = descD->mat;
            double *D = &MAT[tsa(descD, imt, jnt)];
            MAT = descU->mat;
            double *U = &MAT[tsa(descU, imt, jnt)];
            MAT = descV->mat;
            double *V = &MAT[tsa(descV, imt, jnt)];
            MAT = descRk->mat;
            double *Rk = &MAT[tsa(descRk, imt, jnt)];
            int rk = Rk[0];
            printf("%d Tile %d,%d  rk:%d D:%p U:%p V:%p Rk:%p\n", MORSE_My_Mpi_Rank(), imt, jnt, rk, D, U, V, Rk);
            if(imt == jnt){
                _printmat(D, descD->nb, descD->nb, tld(descD));
            }
            else {
                _printmat(U, descD->nb, rk, tld(descU));
                _printmat(V, descD->nb, rk, tld(descV));
            }
        }
    }
}
void _printdescrk(MORSE_desc_t *descZ, int64_t rank){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printf("%d Tile %d,%d  %p\n", MORSE_My_Mpi_Rank(), imt, jnt, A);
            nelm=0;
            for(i=0;i<descZ->nb;i++){
                for(j=0;j<rank;j++){
                    printf("%+.2e ", A[j*tld(descZ)+i]);
                    //printf("%g ", A[j*tld(descZ)+i]);
                    //printf("%g\t", A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
                    nelm++;
                    if(nelm >= _nelm_limit){
                        break;
                    }
                    if(j==_ndim_limit)
                        break;
                }
                if(nelm >= _nelm_limit){
                    printf("\n");
                    break;
                }
                printf("\n");
                if(i==_ndim_limit)
                    break;
            }
        }
    }
}
