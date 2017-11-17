/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * This file contains data structures used in HiCMA.
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 **/

#ifndef __HICMA_STRUCT__
#define __HICMA_STRUCT__

#include "starsh.h"
struct hicma_problem_s;
typedef struct hicma_problem_s HICMA_problem_t;

struct hicma_problem_s {
    STARSH_blrf *starsh_format;

    int ndim;
    double *theta;
    double beta;
    double nu;
    double noise;
};

struct hicma_stat_s;
typedef struct hicma_stat_s HICMA_stat_t;

struct hicma_stat_s {
    int max;
    int min;
    double avg;
};
#endif
