/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file auxdescutil.h
 *
 * This file contains the declarations of auxiliary functions for printing MORSE descriptors..
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#ifndef __AUXDESCUTIL__
#define __AUXDESCUTIL__
#include <stdio.h>
#include "morse.h"
#include "hicma_struct.h"

#define tld(d) (d->mb)
#define tsa(d,i,j) (((j)*(d->mt)+(i))*(d->mb)*(d->nb))

void printmat(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics);
void printmat_format(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics, int format);
void printdescrk(MORSE_desc_t *descZ, int64_t rank);
void printdesc(MORSE_desc_t *descZ);
void _printmat(double * A, int64_t m, int64_t n, int64_t ld);
void _printdescs(MORSE_desc_t *descD,MORSE_desc_t *descU, MORSE_desc_t *descV,  MORSE_desc_t *descRk);
void _printdescrk(MORSE_desc_t *descZ, int64_t rank);

void check_same(MORSE_desc_t *descL, MORSE_desc_t *descR, char diag, char uplo);
void check_same_array(double *L, double *R, int nelm, int line, char *file);

void zget_stat(MORSE_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void zprint_stat(HICMA_stat_t stat);
#endif
