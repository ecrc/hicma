/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * This file contains data structures used in HiCMA.
 *
 * @version 0.1.1
 * @author Kadir Akbudak
 * @date 2018-11-08
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
    double diag;

    // Electrodynamics
    double wave_k;

    int kernel_type;	
    double *point; //coordinates of points
};

struct hicma_stat_s;
typedef struct hicma_stat_s HICMA_stat_t;

struct hicma_stat_s {
    int max;
    int min;
    double avg;
};


struct hicma_context {
    STARSH_blrf *starsh_format;
    char datebuf[128];
    struct tm* tm_info;
    time_t timer;
    int print_progress;   // Print progress about the execution
    int use_fast_hcore_zgemm;
    int store_only_diagonal_tiles;
    int global_check;
    int global_always_fixed_rank;
    int global_fixed_rank;
    int global_omit_computation;
    int num_mpi_ranks;
    int run_potrf;
    int diag_nrows;
    int main_print_index;
    int print_index;
    int print_index_end;
    int main_print_mat;
    int print_mat;
    int use_scratch; // Use scratch memory provided by starpu
    int calc_rank_stat;
};
#endif
