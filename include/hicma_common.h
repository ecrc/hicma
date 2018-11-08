#ifndef __HICMA_COMMON__
#define __HICMA_COMMON__
/**
 * @file hicma_common.h
 *
 * This header file is used inside the library.
 */
#include "starsh.h"
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
#endif
