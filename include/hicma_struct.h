/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * This file contains data structures used in HiCMA.
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/**
 *
 * @file hicma_struct.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon structures
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */

#include <hicma_config.h>
#include <hicma_types.h>
#include <hicma_kernels.h>

#ifndef __HICMA_STRUCT__
#define __HICMA_STRUCT__


#include "starsh.h"

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
    int use_fast_hcore_gemm;
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
typedef struct {
    int id;
    int rank;
} idrank;

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
   
    double reg;
    int numobj;
    int isreg;
    double rad;
    double denst;
    int mesh_points;
    int mordering;
    char* mesh_file;
    char *interpl_file;
    int ntrian;
    int nipp;
};

/**
 * RUNTIME headers to include types of :
 *         - QUARK
 *         - PaRSEC
 *         - StarPU
 */
typedef enum hicma_sched_e {
    HICMA_RUNTIME_SCHED_QUARK,
    HICMA_RUNTIME_SCHED_PARSEC,
    HICMA_RUNTIME_SCHED_STARPU,
} HICMA_sched_t;


/**
 *  Tile matrix descriptor
 *
 *  Matrices are stored in a contiguous data chunk containning in order
 *  A11, A21, A12, A22 with :
 *
 *           n1      n2
 *      +----------+---+
 *      |          |   |    With m1 = lm - (lm%mb)
 *      |          |   |         m2 = lm%mb
 *  m1  |    A11   |A12|         n1 = ln - (ln%nb)
 *      |          |   |         n2 = ln%nb
 *      |          |   |
 *      +----------+---+
 *  m2  |    A21   |A22|
 *      +----------+---+
 *
 */
struct hicma_desc_s;
typedef struct hicma_desc_s HICMA_desc_t;

struct hicma_desc_s {
    // function to get matrix tiles address
    void *(*get_blkaddr)( const HICMA_desc_t*, int, int );
    // function to get matrix tiles leading dimension
    int   (*get_blkldd )( const HICMA_desc_t*, int );
    // function to get matrix tiles MPI rank
    int   (*get_rankof) ( const HICMA_desc_t*, int, int );
    void *mat;        // pointer to the beginning of the matrix
    size_t A21;       // pointer to the beginning of the matrix A21
    size_t A12;       // pointer to the beginning of the matrix A12
    size_t A22;       // pointer to the beginning of the matrix A22
    HICMA_enum styp;  // storage layout of the matrix
    HICMA_enum dtyp;  // precision of the matrix
    int mb;           // number of rows in a tile
    int nb;           // number of columns in a tile
    int bsiz;         // size in elements including padding
    int lm;  	      // number of rows of the entire matrix
    int ln;           // number of columns of the entire matrix
    int lmt;          // number of tile rows of the entire matrix - derived parameter
    int lnt;          // number of tile columns of the entire matrix - derived parameter
    int i;            // row index to the beginning of the submatrix
    int j;            // column index to the beginning of the submatrix
    int m;            // number of rows of the submatrix
    int n;            // number of columns of the submatrix
    int mt;           // number of tile rows of the submatrix - derived parameter
    int nt;           // number of tile columns of the submatrix - derived parameter
    // Data for distributed cases
    int p;            // number of rows of the 2D distribution grid
    int q;            // number of columns of the 2D distribution grid
    int llm;          // number of rows of the 2D distribution grid
    int lln;          // number of columns of the 2D distribution grid
    int llm1;         // number of tile rows of the A11 matrix - derived parameter
    int lln1;         // number of tile columns of the A11 matrix - derived parameter
    int llmt;         // number of tile rows of the local (to a node) matrix
    int llnt;         // number of tile columns of the local (to a node) matrix
    int id;           // identification number of the descriptor
    int occurences;   // identify main matrix desc (occurances=1) or
    // submatrix desc (occurances>1) to avoid unregistering
    // GPU data twice
    int use_mat;      // 1 if we have a pointer to the overall data mat - else 0
    int alloc_mat;    // 1 if we handle the allocation of mat - else 0
    int register_mat; // 1 if we have to register mat - else 0 (handled by the application)
    int myrank;       // MPI rank of the descriptor
    int ooc;          // 1 if the matrix is not to fit in memory
    void *schedopt;   // scheduler (QUARK|StarPU) specific structure
};


/**
 *  HICMA request uniquely identifies each asynchronous function call.
 */
typedef struct hicma_context_s {
    HICMA_sched_t      scheduler;
    int                nworkers;
    int                ncudas;
    int                nthreads_per_worker;
#if defined(HICMA_USE_MPI)
    int                my_mpi_rank;
    int                mpi_comm_size;
#endif
    int                world_size;
    int                group_size;

    /* Boolean flags */
    HICMA_bool         warnings_enabled;
    HICMA_bool         autotuning_enabled;
    HICMA_bool         parallel_enabled;
    HICMA_bool         profiling_enabled;
    HICMA_bool         progress_enabled;

    HICMA_enum         householder;        // "domino" (flat) or tree-based (reduction) Householder
    HICMA_enum         translation;        // In place or Out of place layout conversion

    int                nb;
    int                ib;
    int                nbnbsize;           // tile size in elements (possibly padded)
    int                ibnbsize;           // T or L tile size in elements (---''---)
    int                rhblock;            // block size for tree-based (reduction) Householder
    void              *schedopt;           // structure for runtimes
    int                mpi_outer_init;     // MPI has been initialized outside our functions
} HICMA_context_t;


/**
 *  HICMA request uniquely identifies each asynchronous function call.
 */
typedef struct hicma_request_s {
    HICMA_enum status; // HICMA_SUCCESS or appropriate error code
} HICMA_request_t;


/**
 *  HICMA sequence uniquely identifies a set of asynchronous function calls
 *  sharing common exception handling.
 */
typedef struct hicma_sequence_s {
    HICMA_bool       status;    /* HICMA_SUCCESS or appropriate error code */
    HICMA_request_t *request;   /* failed request                          */
    void            *schedopt;
} HICMA_sequence_t;


/**
 *  HICMA options
 */
typedef struct hicma_option_s {
    HICMA_sequence_t *sequence;
    HICMA_request_t  *request;
    int               profiling;
    int               parallel;
    int               priority;
    int               nb;
    size_t            ws_wsize;
    size_t            ws_hsize;
    void             *ws_worker;  /*> Workspace located on the worker        */
    void             *ws_host;    /*> Workspace *always* located on the host */
    void             *schedopt;
} HICMA_option_t;

#endif
