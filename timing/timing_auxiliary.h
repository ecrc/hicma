/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file hicma_d.h
 *
 *  HiCMA computational routines
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

#include "hicma_struct.h"

#ifdef __cplusplus
extern "C" {
#endif
#define PROBLEM_TYPE_RND    1
#define PROBLEM_TYPE_SS     2
#define PROBLEM_TYPE_RNDUSR 3
#define PROBLEM_TYPE_FILE   4
#define PROBLEM_TYPE_GEOSTAT   5
#define PROBLEM_TYPE_EDSIN  6
#define PROBLEM_TYPE_GEOSTAT_POINT   7
#define PROBLEM_TYPE_ST_3D_EXP  8
#define PROBLEM_TYPE_ST_3D_SQEXP  9
#define PROBLEM_TYPE_3D_RBF_VIRUS  12
#define PROBLEM_TYPE_3D_RBF_CUBE  13
#define PROBLEM_TYPE_AC_3D  14
#define PROBLEM_TYPE_ST_2D_EXP  15
#define PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS_BIVARIATE 108
#define PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT 109
#define PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS2_BIVARIATE 110
#define PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT 111

void generate_problem(
        int probtype, //problem type defined in hicma_constants.h 
        char sym,     // symmetricity of problem: 'N' or 'S'
        double decay, // decay of singular values. Will be used in PROBLEM_TYPE_RND. Set 0 for now.
        int _M,       // number of rows/columns of matrix
        int _nb,      // number of rows/columns of a single tile
        int _mt,      // number of tiles in row dimension
        int _nt,      // number of tiles in column dimension
        HICMA_problem_t *hicma_problem // pointer to hicma struct (starsh format will be used to pass coordinate info to number generation and compression phase)
        );
#ifdef __cplusplus
}
#endif
