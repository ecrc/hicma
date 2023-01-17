/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file dproblem.c
 *
 * This file contains the function for generating a problem.
 * This problem can then be used to generate exact (dense) or
 * approximate matrix representing the problem.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2019-11-21
 **/

#include "hicma.h"
#include "starsh.h"
#include "starsh-spatial.h"
#include "starsh-electrodynamics.h"
#include "starsh-randtlr.h"
#include "starsh-rbf.h"
#include "starsh-acoustic.h"

#include "timing_auxiliary.h"

#include <stdlib.h>

void generate_problem(int probtype, char sym, double decay, int _M, int block_size, int _mt, int _nt,
                      HICMA_problem_t *hicma_problem) {
    char dtype = 'd';
    STARSH_kernel *kernel;
    STARSH_int shape[] = {_M, _M};
    if (_mt != _nt) {
        fprintf(stderr, "number of tiles must be equal in each dimension: mt:%d nt:%d\n", _mt, _nt);
    }
    int info;
    STARSH_problem *problem = NULL;
    void *data;
    int ndim = 2; // two dimensional matrix is used all time
    char *strproblem;
    if (probtype == PROBLEM_TYPE_RND) {
        int kernel_type = STARSH_RANDTLR_KERNEL1;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_RANDTLR, kernel_type, STARSH_RANDTLR_NB, block_size,
                                  STARSH_RANDTLR_DECAY, decay, STARSH_RANDTLR_DIAG, hicma_problem->noise,
                                  0);
        strproblem = "Randomly generated matrix";
    } else if (probtype == PROBLEM_TYPE_SS) {
        int kernel_type = STARSH_SPATIAL_SQREXP_SIMD;
        srand(0); // FIXME
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->beta, STARSH_SPATIAL_NU, hicma_problem->nu,
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, place,
                                  0);
        strproblem = "ST_2D_SQEXP";
    } else if (probtype == PROBLEM_TYPE_ST_2D_EXP) {
        int kernel_type = STARSH_SPATIAL_EXP_SIMD;
        srand(0); // FIXME
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->beta, STARSH_SPATIAL_NU, hicma_problem->nu,
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, place,
                                  0);
        strproblem = "ST_2D_EXP";
    } else if (probtype == PROBLEM_TYPE_ST_3D_EXP) {
        int kernel_type = STARSH_SPATIAL_EXP_SIMD;
        srand(0); // FIXME
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        int probndim = 3;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, probndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->beta, STARSH_SPATIAL_NU, hicma_problem->nu,
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, place,
                                  0);
        strproblem = "ST_3D_EXP";
    } else if (probtype == PROBLEM_TYPE_ST_3D_SQEXP) {
        int kernel_type = STARSH_SPATIAL_SQREXP_SIMD;
        srand(0); // FIXME
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        int probndim = 3;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, probndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->beta, STARSH_SPATIAL_NU, hicma_problem->nu,
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, place,
                                  0);
        strproblem = "ST_3D_SQEXP";
    } else if (probtype == PROBLEM_TYPE_EDSIN) {
        int kernel_type = STARSH_ELECTRODYNAMICS_SIN;
        srand(0); // FIXME
        enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
        double wave_k = hicma_problem->wave_k;
        double diag = hicma_problem->diag;
        if (diag == 0) {
            info = starsh_application((void **) &data, &kernel, _M, dtype,
                                      STARSH_ELECTRODYNAMICS, kernel_type, STARSH_ELECTRODYNAMICS_NDIM, ndim,
                                      STARSH_ELECTRODYNAMICS_K, wave_k,
                                      STARSH_ELECTRODYNAMICS_PLACE, place,
                                      0);
        } else {
            info = starsh_application((void **) &data, &kernel, _M, dtype,
                                      STARSH_ELECTRODYNAMICS, kernel_type, STARSH_ELECTRODYNAMICS_NDIM, ndim,
                                      STARSH_ELECTRODYNAMICS_K, wave_k,
                                      STARSH_ELECTRODYNAMICS_DIAG, diag,
                                      STARSH_ELECTRODYNAMICS_PLACE, place,
                                      0);
        }
        strproblem = "ELECTRODYNAMICS_SINUS";
    } else if (probtype == PROBLEM_TYPE_GEOSTAT) {
        srand(0);
        int kernel_type = hicma_problem->kernel_type;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->theta[1], STARSH_SPATIAL_NU,
                                  hicma_problem->theta[2],
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE1, //1,
                                  STARSH_SPATIAL_SIGMA, hicma_problem->theta[0],
                                  0);
        strproblem = "GEOSTAT";
    } else if (probtype == PROBLEM_TYPE_GEOSTAT_POINT) {
        int kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
        info = starsh_ssdata_init((STARSH_ssdata **) &data, _M, ndim, hicma_problem->point, hicma_problem->theta[1],
                                  hicma_problem->theta[2], hicma_problem->noise, hicma_problem->theta[0]);
        strproblem = "GEOSTAT_POINT";
    } else if (probtype == PROBLEM_TYPE_3D_RBF_VIRUS) {
        int kernel_type = hicma_problem->kernel_type;
        double reg = hicma_problem->reg;
        int isreg = hicma_problem->isreg;
        double rad = hicma_problem->rad;
        double denst = hicma_problem->denst;
        int mesh_points = hicma_problem->mesh_points;
        int ordering = hicma_problem->mordering;
        char *mesh_file = hicma_problem->mesh_file;
        int numobj = hicma_problem->numobj;
        int problem_ndim = 3;
        info = starsh_generate_3d_rbf_mesh_coordinates_virus((STARSH_mddata **) &data, mesh_file, mesh_points,
                                                             problem_ndim, kernel_type, numobj, isreg, reg, rad, denst,
                                                             ordering);
        kernel = starsh_generate_3d_virus;
        strproblem = "SARS-CoV-2";
    } else if (probtype == PROBLEM_TYPE_3D_RBF_CUBE) {
        int kernel_type = hicma_problem->kernel_type;
        double reg = hicma_problem->reg;
        int isreg = hicma_problem->isreg;
        double rad = hicma_problem->rad;
        int mesh_points = hicma_problem->mesh_points;
        int ordering = hicma_problem->mordering;
        char *mesh_file = hicma_problem->mesh_file;
        int problem_ndim = 3;
        info = starsh_generate_3d_rbf_mesh_coordinates_cube((STARSH_mddata **) &data, mesh_points, problem_ndim,
                                                            kernel_type, isreg, reg, rad, ordering);
        kernel = starsh_generate_3d_cube;
        strproblem = "Cube";
    } else if (probtype == PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS_BIVARIATE) {
        srand(0);
        int kernel_type = hicma_problem->kernel_type;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->theta[2], STARSH_SPATIAL_NU,
                                  hicma_problem->theta[3],
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE3, //1,
                                  STARSH_SPATIAL_SIGMA, hicma_problem->theta[0],
                                  STARSH_SPATIAL_SIGMA2, hicma_problem->theta[1],
                                  STARSH_SPATIAL_NU2, hicma_problem->theta[4],
                                  STARSH_SPATIAL_CORR, hicma_problem->theta[5],
                                  0);
        strproblem = "GEOSTAT_PARSIMONIOUS_BIVARIATE";
    } else if (probtype == PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS2_BIVARIATE) {
        srand(0);
        int kernel_type = hicma_problem->kernel_type;
        info = starsh_application((void **) &data, &kernel, _M, dtype,
                                  STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                  STARSH_SPATIAL_BETA, hicma_problem->theta[2], STARSH_SPATIAL_NU,
                                  hicma_problem->theta[3],
                                  STARSH_SPATIAL_NOISE, hicma_problem->noise,
                                  STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE4, //1,
                                  STARSH_SPATIAL_SIGMA, hicma_problem->theta[0],
                                  STARSH_SPATIAL_SIGMA2, hicma_problem->theta[1],
                                  STARSH_SPATIAL_NU2, hicma_problem->theta[4],
                                  STARSH_SPATIAL_CORR, hicma_problem->theta[5],
                                  0);
        strproblem = "GEOSTAT_PARSIMONIOUS2_BIVARIATE";
    } else if (probtype == PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT) {
        int kernel_type = STARSH_SPATIAL_PARSIMONIOUS_SIMD;
        info = starsh_ssdata_init_parsimonious((STARSH_ssdata **) &data, _M, ndim, hicma_problem->point,
                                               hicma_problem->theta[0],
                                               hicma_problem->theta[1], hicma_problem->theta[2],
                                               hicma_problem->theta[3],
                                               hicma_problem->theta[4], hicma_problem->theta[5], hicma_problem->noise);
        strproblem = "GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT";
        if (info != 0) {
            printf("wrong parameters for starsh_ssdata_init\n");
        }
        info = starsh_ssdata_get_kernel(&kernel, data, kernel_type);
        if (info != 0) {
            printf("wrong parameters for starsh_ssdata_get_kernel\n");
        }
    } else if (probtype == PROBLEM_TYPE_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT) {
        int kernel_type = STARSH_SPATIAL_PARSIMONIOUS2_SIMD;
        info = starsh_ssdata_init_parsimonious((STARSH_ssdata **) &data, _M, ndim, hicma_problem->point,
                                               hicma_problem->theta[0],
                                               hicma_problem->theta[1], hicma_problem->theta[2],
                                               hicma_problem->theta[3],
                                               hicma_problem->theta[4], hicma_problem->theta[5], hicma_problem->noise);
        strproblem = "GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT";
        if (info != 0) {
            printf("wrong parameters for starsh_ssdata_init\n");
        }
        info = starsh_ssdata_get_kernel(&kernel, data, kernel_type);
        if (info != 0) {
            printf("wrong parameters for starsh_ssdata_get_kernel\n");
        }
    } else if (probtype == PROBLEM_TYPE_AC_3D) {
        dtype = 'z';
        int trian = hicma_problem->ntrian;
        int nipp = hicma_problem->nipp;
        int mesh_points = hicma_problem->mesh_points;
        int ordering = hicma_problem->mordering;
        int problem_ndim = 3;
        char *mesh_file = hicma_problem->mesh_file;
        char *interpl_file = hicma_problem->interpl_file;

        info = starsh_generate_3d_acoustic_coordinates((STARSH_acdata **) &data, mesh_points, problem_ndim, trian, nipp,
                                                       0, mesh_file, interpl_file);
        kernel = starsh_generate_3d_acoustic;
        strproblem = "Acoustic Scattering";
    } else {
        fprintf(stderr, "Unknown type of STARS-H problem:%d. Exiting...\n", probtype);
    }
    if (info != 0) {
        printf("wrong parameters for starsh_application()\n");
        exit(info);
    }
    starsh_problem_new(&problem, ndim, shape, sym, dtype, data, data, kernel, strproblem);
    STARSH_cluster *cluster;
    info = starsh_cluster_new_plain(&cluster, data, _M, block_size);
    if (info != 0) {
        printf("Error in creation of cluster\n");
        exit(info);
    }
    STARSH_blrf *F;
    info = starsh_blrf_new_tlr(&F, problem, sym, cluster, cluster);
    if (info != 0) {
        printf("Error in creation of format\n");
        exit(info);
    }
    hicma_problem->starsh_format = F;
    HICMA_set_starsh_format(hicma_problem->starsh_format);
}
