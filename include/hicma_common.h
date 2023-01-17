/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

#ifndef __HICMA_COMMON__
#define __HICMA_COMMON__
/**
 * @file hicma_common.h
 *
 * This header file is used inside the library.
 */
#include <hicma.h>
#include <hicma_struct.h>
#include <hicma_z.h>
#include <control/hicma_compute_z.h>
#include <control/hicma_compute_s.h>
#include <control/hicma_compute_d.h>
#include <control/hicma_compute_c.h>

#define PROGRESS(str) \
    if(print_progress){ \
        int myrank = HICMA_My_Mpi_Rank();\
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S",tm_info); \
        fprintf(stderr, "%d:%s\t%d\t%s\t%s\n", myrank, datebuf, __LINE__, __func__, str);\
        fflush(stderr);\
    }
//#undef PROGRESS
//#define PROGRESS(str)
extern struct hicma_context hicma_context;

void printdescrk(HICMA_desc_t *descZ, int64_t rank);
void printdesc(HICMA_desc_t *descZ);
void _printdescs(HICMA_desc_t *descD,HICMA_desc_t *descU, HICMA_desc_t *descV,  HICMA_desc_t *descRk);
void _printdescrk(HICMA_desc_t *descZ, int64_t rank);

void check_same(HICMA_desc_t *descL, HICMA_desc_t *descR, char diag, char uplo);
void dget_stat(HICMA_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void dprint_stat(HICMA_stat_t stat);
void zget_stat(HICMA_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void zprint_stat(HICMA_stat_t stat);
int HICMA_Lapack_to_Tile(void *Af77, int LDA, HICMA_desc_t *A);
int HICMA_Tile_to_Lapack(HICMA_desc_t *A, void *Af77, int LDA);
#endif
