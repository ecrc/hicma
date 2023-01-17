/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <hicma_common.h>
#include <misc/auxdescutil.h>

int64_t hicma_nelm_limit = 100000;
int64_t hicma_nrow_limit = 100;
int64_t hicma_ncol_limit = 100;
#define hicma_tld(d) (d->mb)
#define hicma_tsa(d,i,j) (((j)*(d->mt)+(i))*(d->mb)*(d->nb))
int64_t hicma_2_hicma_nelm_limit = 1200;
int64_t hicma_2_nrows_limit = 37;
int64_t hicma_2_ncols_limit = 33;
int64_t hicma_2_ndim_limit = 37;

void printdescrk(HICMA_desc_t *descZ, int64_t rank){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            printf("Tile %d,%d\n", imt, jnt);
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[hicma_tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printmat(A, descZ->nb, rank, hicma_tld(descZ), -1, -1);
        }
    }
}
void printdesc(HICMA_desc_t *descZ){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            printf("Tile %d,%d\n", imt, jnt);
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[hicma_tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printmat(A, descZ->nb, descZ->nb, hicma_tld(descZ), -1, -1);
        }
    }
}
void _printdescs(HICMA_desc_t *descD,HICMA_desc_t *descU, HICMA_desc_t *descV,  HICMA_desc_t *descRk){
    int64_t i, j, imt, jnt;
    printf("\n");
    for(imt=0;imt<descD->mt;imt++){
        for(jnt=0;jnt<descD->nt;jnt++){
            if(imt < jnt) continue;
            double *MAT = descD->mat;
            double *D = &MAT[hicma_tsa(descD, imt, jnt)];
            MAT = descU->mat;
            double *U = &MAT[hicma_tsa(descU, imt, jnt)];
            MAT = descV->mat;
            double *V = &MAT[hicma_tsa(descV, imt, jnt)];
            MAT = descRk->mat;
            double *Rk = &MAT[hicma_tsa(descRk, imt, jnt)];
            int rk = Rk[0];
            printf("%d Tile %d,%d  rk:%d D:%p U:%p V:%p Rk:%p\n", HICMA_My_Mpi_Rank(), imt, jnt, rk, D, U, V, Rk);
            if(imt == jnt){
                _printmat(D, descD->nb, descD->nb, hicma_tld(descD));
            }
            else {
                _printmat(U, descD->nb, rk, hicma_tld(descU));
                _printmat(V, descD->nb, rk, hicma_tld(descV));
            }
        }
    }
}
void _printdescrk(HICMA_desc_t *descZ, int64_t rank){
    double *MAT = descZ->mat;
    int64_t i, j, imt, jnt, nelm = 0;
    printf("\n");
    for(imt=0;imt<descZ->mt;imt++){
        for(jnt=0;jnt<descZ->nt;jnt++){
            //double *A = &MAT[(jnt*descZ->mt+imt)*descZ->mb*descZ->nb];
            double *A = &MAT[hicma_tsa(descZ, imt, jnt)];
            //double *A =(double*) RTBLKADDR(descZ, double, imt, jnt); // does not work
            printf("%d Tile %d,%d  %p\n", HICMA_My_Mpi_Rank(), imt, jnt, A);
            nelm=0;
            for(i=0;i<descZ->nb;i++){
                for(j=0;j<rank;j++){
                    printf("%+.2e ", A[j*hicma_tld(descZ)+i]);
                    //printf("%g ", A[j*hicma_tld(descZ)+i]);
                    //printf("%g\t", A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->mb,A[j*descZ->mb+i]);
                    //printf("(%d,%d,%d) %g\t", i,j,descZ->n,A[j*descZ->n+i]);
                    //printf("(%d,%d,%d) %g [%d %d %d]\t", i,j,descZ->n,A[j*descZ->n+i], descZ->m, descZ->lm, descZ->ln);
                    nelm++;
                    if(nelm >= hicma_2_hicma_nelm_limit){
                        break;
                    }
                    if(j==hicma_2_ndim_limit)
                        break;
                }
                if(nelm >= hicma_2_hicma_nelm_limit){
                    printf("\n");
                    break;
                }
                printf("\n");
                if(i==hicma_2_ndim_limit)
                    break;
            }
        }
    }
}
void check_same(HICMA_desc_t *descL, HICMA_desc_t *descR, char diag, char uplo){
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
