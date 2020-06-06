#ifndef __HICMA_COMMON__
#define __HICMA_COMMON__
/**
 * @file hicma_common.h
 *
 * This header file is used inside the library.
 */
#include "morse.h"
#include "hicma_struct.h"
#define PROGRESS(str) \
    if(print_progress){ \
        int myrank = MORSE_My_Mpi_Rank();\
        time(&timer); \
        tm_info = localtime(&timer); \
        strftime(datebuf, 26, "%Y-%m-%d %H:%M:%S",tm_info); \
        fprintf(stderr, "%d:%s\t%d\t%s\t%s\n", myrank, datebuf, __LINE__, __func__, str);\
        fflush(stderr);\
    }
//#undef PROGRESS
//#define PROGRESS(str)
extern struct hicma_context hicma_context;

void printdescrk(MORSE_desc_t *descZ, int64_t rank);
void printdesc(MORSE_desc_t *descZ);
void _printdescs(MORSE_desc_t *descD,MORSE_desc_t *descU, MORSE_desc_t *descV,  MORSE_desc_t *descRk);
void _printdescrk(MORSE_desc_t *descZ, int64_t rank);

void check_same(MORSE_desc_t *descL, MORSE_desc_t *descR, char diag, char uplo);
void zget_stat(MORSE_enum uplo, double *Ark, size_t m, size_t n, size_t ld,  HICMA_stat_t *stat);
void zprint_stat(HICMA_stat_t stat);
#endif
