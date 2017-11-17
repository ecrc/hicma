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
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/
#ifndef __AUXDESCUTIL__
#define __AUXDESCUTIL__
#include <stdio.h>
#include "morse.h"
#include "hicma_struct.h"

#define tld(d) (d->mb)
#define tsa(d,i,j) (((j)*(d->mt)+(i))*(d->mb)*(d->nb))

void printmat(double * A, int64_t m, int64_t n, int64_t ld, int irs, int ics);
void printdescrk(MORSE_desc_t *descZ, int64_t rank);
void printdesc(MORSE_desc_t *descZ);

void zget_stat(MORSE_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void zprint_stat(HICMA_stat_t stat);
#endif
