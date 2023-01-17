/**
 * @copyright (c) 2012 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file zproblem.c
 *
 * This file contains the function for generating a problem.
 * This problem can then be used to generate exact (dense) or
 * approximate matrix representing the problem.
 *
 * HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/
#include <math.h>
#include <sys/time.h>
#include "starsh.h"
#include "misc/auxcompute_z.h"
#include "starsh-spatial.h"
#include "starsh-electrodynamics.h"
#include "starsh-randtlr.h"
#include "hicma_constants.h"
#include "hicma_struct.h"
#include "hicma_common.h"
#include <assert.h>
#include "misc/auxdescutil.h"
#include "hicma.h"

int print_starsh_info = 1;

void __generate(
        int probtype,
        char sym,
        double decay,
        int _M,
        int _nb,
        int _mt,
        int _nt,
        int num_descs,
        HICMA_desc_t **descs,
        HICMA_desc_t **descsU,
        HICMA_desc_t **descsV,
        HICMA_desc_t **descsD,
        HICMA_desc_t **descsrk,
        double **lrarrays,
        int *pactual_maxrank, double *pactual_avgrank,
        int diag_dense,
        STARSH_blrf **mpiF,
        double *initial_theta,
        int ndim, double beta, double nu, double noise,
        int kernel_type, double *point,
        HICMA_problem_t *hicma_problem
) {
    int idesc;
    int print_prep_time = 0; /// 1:print timing info preprocessing, 0: no print
    int starsh_printmatelm = 0;
    int starsh_printmatelmD = 0;
    int starsh_printmatelmU = 0;
    int starsh_printmatelmV = 0;
    struct timeval t1s, t1e;
    gettimeofday(&t1s, 0);
    for (idesc = 0; idesc < num_descs; idesc++) { //foreachMatrix
        //HICMA_desc_t *descX   = descs  [idesc];
        HICMA_desc_t *descXU = descsU[idesc];
        HICMA_desc_t *descXV = descsV[idesc];
        HICMA_desc_t *descXD = descsD[idesc];
        HICMA_desc_t *descXrk = descsrk[idesc];
        double *lrarray_double = NULL; //lrarrays[idesc];

        char calcAcc = 'c';            /// (c)alculate accuracy or (s)kip the calculation
        int seed = 0;
        double *array_double = NULL;
        int starsh_printmat = 1;
        int stars_print = 0;
        int stars_free = 1;
        int stars_check_input = 0;
        int sqrtm = sqrt(_M);
        int block_size = _nb;
        char dtype = 'd';
        //STARSH_SPATIAL_MATERN2_SIMD; //sameh
        STARSH_kernel *kernel;
        //STARSH_cluster *C = NULL;
        //STARSH_problem *P = NULL;
        STARSH_int shape[] = {_M, _M};
        if (_mt != _nt) {
            fprintf(stderr, "number of tiles must be equal in each dimension: mt:%d nt:%d\n", _mt, _nt);
        }
        //assert(_mt == _nt);   //check if the number of tiles in both dims are same
        //assert((_M/_nb) == _mt); //check if par _mt is equal to computed #tiles


        double sigma = 1.0; //FIXME TODO KADIR ADD !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        if (probtype == HICMA_STARSH_PROB_RND) {
            int info;
            // Since there is only one kernel for rndtiled, kernel_type is ignored
            int kernel_type = STARSH_RANDTLR_KERNEL1;
            STARSH_cluster *C = NULL;
            STARSH_problem *P = NULL;

            if (print_starsh_info) {
                printf("M:%d block_size:%d decay:%g diag:%g\n", _M, block_size, decay, noise);
            }

            STARSH_randtlr *data;
            // STARSH_RANDTLR for random tile low-rank matrix
            // STARSH_RANDTLR_NB to indicate next parameter shows size of tile
            // STARSH_RANDTLR_DECAY to indicate next parameter is decay of singular
            //   values
            // STARSH_RANDTLR_DIAG to indicate next parameter is addition for diagonal
            //   elements
            // 0 at the end to indicate end of arguments
            info = starsh_application((void **) &data, &kernel, _M, dtype,
                                      STARSH_RANDTLR, kernel_type, STARSH_RANDTLR_NB, block_size,
                                      STARSH_RANDTLR_DECAY, decay, STARSH_RANDTLR_DIAG, noise,
                                      0);
            if (info != 0) {
                printf("wrong parameters for random tile low-rank matrix\n");
                exit(info);
            }


            // Init problem with given data and kernel and print short info
            starsh_problem_new(&P, 2, shape, sym, 'd', data, data, kernel, "Randomly generated matrix");

            if (print_starsh_info) {
                printf("\nDecay:%e\n", decay);
                starsh_problem_info(P);
            }
            // Init tiled cluster for tiled low-rank approximation and print info
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, data, _M, _nb);
            if (info != 0) {
                printf("Error in creation of cluster\n");
                exit(info);
            }
            if (print_starsh_info) {
                starsh_cluster_info(cluster);
            }
            STARSH_blrf *F;
            info = starsh_blrf_new_tlr(&F, P, sym, cluster, cluster);
            if (info != 0) {
                printf("Error in creation of format\n");
                exit(info);
            }
            if (print_starsh_info) {
                starsh_blrf_info(F);
            }
            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_SS) {
            fprintf(stderr, "%s %d Spatial Statistics Application with SqExp Kernel: beta:%g nu:%g noise:%g\n",
                    __FILE__, __LINE__,
                    beta,
                    nu,
                    noise);
            // Correlation length
            // double beta;
            // Smoothing parameter for Matern kernel
            //double nu = 0.5;
            // Set level of noise
            //double noise = 1.e-4;

            //for small blocks size, for testing, low rank
            //double noise = _M  1 ... ;
            // Possible values for kernel_type are:
            //      STARSH_SPATIAL_EXP, STARSH_SPATIAL_EXP_SIMD
            //      STARSH_SPATIAL_SQREXP, STARSH_SPATIAL_SQREXP_SIMD
            //      STARSH_SPATIAL_MATERN, STARSH_SPATIAL_MATERN_SIMD
            //      STARSH_SPATIAL_MATERN2, STARSH_SPATIAL_MATERN2_SIMD
            int kernel_type = STARSH_SPATIAL_SQREXP_SIMD;
            srand(0); // FIXME
            /*if((sqrtm*sqrtm) != _M){*/
            /*fprintf(stderr,"M must be square of integer.\n");*/
            /*exit(1);*/
            /*}*/
            /*printf("\n\nSPATIAL STATISTICS\n\n");*/
            //fit a uniform grid on the given size
            enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
            //you have grid and move particles
            //enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_OBSOLETE1;
            STARSH_ssdata *data;
            int info;
            // STARSH_SPATIAL for spatial statistics problem
            // kernel_type is enum type, for possible values llok into starsh-spatial.h
            // STARSH_SPATIAL_NDIM to indicate next parameter shows dimensionality of
            //   spatial statistics problem
            // STARSH_SPATIAL_BETA to indicate next parameter is correlation length
            // STARSH_SPATIAL_NU to indicate next parameter is smoothing parameter for
            //   Matern kernel
            // STARSH_SPATIAL_NOISE to indicate next parameter is a noise
            // 0 at the end to indicate end of arguments
            info = starsh_application((void **) &data, &kernel, _M, dtype,
                                      STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                      STARSH_SPATIAL_BETA, beta, STARSH_SPATIAL_NU, nu,
                                      STARSH_SPATIAL_NOISE, noise,
                                      STARSH_SPATIAL_PLACE, place,
                                      0);
            if (info != 0) {
                printf("wrong parameters for spatial statistics problem\n");
                return;
            }
            int ndim = 2;
            // Init problem with given data and kernel and print short info
            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, sym, dtype, data, data,
                                      kernel, "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            if (print_starsh_info) {
                starsh_problem_info(problem);
            }
            //printf("STARSH problem was succesfully generated\n");
            //starsh_problem_info(problem);
            // Set up clusterization (divide rows and columns into blocks)
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, data, _M, block_size);
            if (info != 0) {
                printf("Error in creation of cluster\n");
                exit(info);
            }
            if (print_starsh_info) {
                starsh_cluster_info(cluster);
            }
            STARSH_blrf *F;
            info = starsh_blrf_new_tlr(&F, problem, sym, cluster, cluster);
            if (info != 0) {
                printf("Error in creation of format\n");
                exit(info);
            }
            if (print_starsh_info) {
                starsh_blrf_info(F);
            }

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_EDSIN) {
            int kernel_type = STARSH_ELECTRODYNAMICS_SIN;
            srand(0); // FIXME
            enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_UNIFORM;
            //you have grid and move particles
            //enum STARSH_PARTICLES_PLACEMENT place = STARSH_PARTICLES_OBSOLETE1;
            // Wave number, >= 0
            double wave_k = hicma_problem->wave_k;
            double diag = hicma_problem->diag;
            fprintf(stderr, "%s %d Electro Dynamics Application with Sinus Kernel wave_k:%g diag:%g\n", __FILE__,
                    __LINE__,
                    hicma_problem->wave_k,
                    hicma_problem->diag
            );
            STARSH_eddata *data;
            int info;
            // STARSH_ELECTRODYNAMICS for electrodynamics problem
            // kernel_type is enum type, for possible values look into starsh-electrodynamics.h
            // STARSH_ELECTRODYNAMICS_NDIM to indicate next parameter shows dimensionality of
            //   spatial statistics problem
            // STARSH_ELECTRODYNAMICS_K to indicate next parameter is wave number
            // STARSH_ELECTRODYNAMICS_DIAG to indicate next parameter is diagonal values
            // 0 at the end to indicate end of arguments
            if (print_starsh_info) {
                printf("\n_M:%d wave_k=%g\n", _M, wave_k);
            }
            info = starsh_application((void **) &data, &kernel, _M, dtype,
                                      STARSH_ELECTRODYNAMICS, kernel_type, STARSH_ELECTRODYNAMICS_NDIM, ndim,
                                      STARSH_ELECTRODYNAMICS_K, wave_k, STARSH_ELECTRODYNAMICS_DIAG, diag,
                                      STARSH_ELECTRODYNAMICS_PLACE, place,
                                      0);
            if (info != 0) {
                printf("wrong parameters for electrodynamics problem\n");
                return;
            }
            int ndim = 2;
            // Init problem with given data and kernel and print short info
            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, sym, dtype, data, data,
                                      kernel, "Electrodynamics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            if (print_starsh_info) {
                printf("\nDecay:%e\n", decay);
                starsh_problem_info(problem);
            }
            //printf("STARSH problem was succesfully generated\n");
            //starsh_problem_info(problem);
            // Set up clusterization (divide rows and columns into blocks)
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

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_RNDUSR) {
            assert(0 == "Not ready");
            int nblocks = _mt;
            int blocksize = _nb;
            //double noise = 0.0;
            //double noise = _acc/decay;
            //printf("decay:%.2e _acc:%.2e noise:%.2e\n", decay, _acc, noise);
            int n = _M;
            STARSH_randtlr *data;
            //starsh_rndtiled_gen(&data, &kernel, nblocks, block_size, decay, noise);
            // Init problem with given data and kernel and print short info
            STARSH_problem *P;
            starsh_problem_new(&P, 2, shape, sym, 'd', data, data, kernel,
                               "Randomly generated matrix");
            /*starsh_problem_info(P);*/
            // Create new problem out of dense matrix
            Array *A;
            starsh_problem_to_array(P, &A);
            double *matrix = A->data;
            int i;
            for (i = 0; i < n; i++)
                matrix[i * (n + 1)] += 1;
            /*printmat(matrix, _M, _M, _M);*/
            starsh_problem_free(P);
            starsh_problem_from_array(&P, A, sym);
            // Init tiled cluster for tiled low-rank approximation and print info
            STARSH_cluster *C;
            //starsh_cluster_new_plain(&C, data, n, block_size);
            starsh_cluster_new_plain(&C, A, n, block_size);
            /*starsh_cluster_info(C);*/
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT) {
            /*fprintf(stderr, "%s %d Geostat Application sigma:%g beta:%g nu:%g noise:%g\n", __FILE__, __LINE__, 
                    initial_theta[0],
                    initial_theta[1],
                    initial_theta[2],
                    noise
                    );
	    */
            //double initial_theta[3] = {1, 0.1, 0.5};
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            srand(0);
            //int kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
            //double noise = 0;
            STARSH_ssdata *ssdata;
            info = starsh_application((void **) &ssdata, &kernel, _M, dtype,
                                      STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                      STARSH_SPATIAL_BETA, initial_theta[1], STARSH_SPATIAL_NU, initial_theta[2],
                                      STARSH_SPATIAL_NOISE, noise,
                                      STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE1, //1,
                                      STARSH_SPATIAL_SIGMA, initial_theta[0],
                                      0);
            if (info != 0) {
                printf("wrong parameters for spatial statistics problem\n");
            }


            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            //    printf("STARSH problem was succesfully generated\n");
            //    starsh_problem_info(problem);
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_POINT) {
            //   printf("Geostat");
            //double initial_theta[3] = {1, 0.1, 0.5};
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            //srand(0);  //no need in case og real dataset
            int kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
            //double noise = 0;
            STARSH_ssdata *ssdata;
//int starsh_ssdata_init(STARSH_ssdata **data, STARSH_int count, int ndim,
            /*double *point, double beta, double nu, double noise, double sigma)*/
            info = starsh_ssdata_init(&ssdata, _M, ndim, point, initial_theta[1],
                                      initial_theta[2], noise, initial_theta[0]);
            /*info = starsh_application((void **)&ssdata, &kernel, _M, dtype,*/
            /*STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,*/
            /*STARSH_SPATIAL_BETA, initial_theta[1], STARSH_SPATIAL_NU, initial_theta[2],*/
            /*STARSH_SPATIAL_NOISE, noise,*/
            /*STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE1, //1,*/
            /*STARSH_SPATIAL_SIGMA, initial_theta[0],*/
            /*0);*/
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_init\n");
            }
            /*int starsh_ssdata_get_kernel(STARSH_kernel **kernel, STARSH_ssdata *data,*/
            /*enum STARSH_SPATIAL_KERNEL type)*/
            info = starsh_ssdata_get_kernel(&kernel, ssdata, kernel_type);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_get_kernel\n");
            }


            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            //    printf("STARSH problem was succesfully generated\n");
            //    starsh_problem_info(problem);
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE) {
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            srand(0);
            STARSH_ssdata *ssdata;
            info = starsh_application((void **) &ssdata, &kernel, _M, dtype,
                                      STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                      STARSH_SPATIAL_BETA, initial_theta[2], STARSH_SPATIAL_NU, initial_theta[3],
                                      STARSH_SPATIAL_NOISE, noise,
                                      STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE3, //1,
                                      STARSH_SPATIAL_SIGMA, initial_theta[0],
                                      STARSH_SPATIAL_SIGMA2, initial_theta[1],
                                      STARSH_SPATIAL_NU2, initial_theta[4],
                                      STARSH_SPATIAL_CORR, initial_theta[5],
                                      0);
            if (info != 0)
                printf("wrong parameters for spatial statistics problem\n");
            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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
            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE) {
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            srand(0);

            STARSH_ssdata *ssdata;
            info = starsh_application((void **) &ssdata, &kernel, _M, dtype,
                                      STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                      STARSH_SPATIAL_BETA, initial_theta[2], STARSH_SPATIAL_NU, initial_theta[3],
                                      STARSH_SPATIAL_NOISE, noise,
                                      STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE4, //1,
                                      STARSH_SPATIAL_SIGMA, initial_theta[0],
                                      STARSH_SPATIAL_SIGMA2, initial_theta[1],
                                      STARSH_SPATIAL_NU2, initial_theta[4],
                                      STARSH_SPATIAL_CORR, initial_theta[5],
                                      0);
            if (info != 0)
                printf("wrong parameters for spatial statistics problem\n");

            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");

            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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
            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN) {
            //double initial_theta[3] = {1, 0.1, 0.5};
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            srand(0);
            //int kernel_type = STARSH_SPATIAL_MATERN2_SIMD;
            int kernel_type = STARSH_SPATIAL_NON_GAUSSIAN_SIMD;
            double noise = 0;

            STARSH_ssdata *ssdata;
            info = starsh_application((void **) &ssdata, &kernel, _M, dtype,
                                      STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,
                                      STARSH_SPATIAL_BETA, initial_theta[0], STARSH_SPATIAL_NU, initial_theta[1],
                                      STARSH_SPATIAL_NOISE, noise,
                                      STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE1, //1,
                                      STARSH_SPATIAL_SIGMA, 1.0,
                                      0);

            if (info != 0) {
                printf("wrong parameters for spatial statistics problem\n");
            }


            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            //    printf("STARSH problem was succesfully generated\n");
            //    starsh_problem_info(problem);
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_NON_GAUSSIAN_POINT) {
            //   printf("Geostat");
            //double initial_theta[3] = {1, 0.1, 0.5};
            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            //srand(0);  //no need in case og real dataset
            int kernel_type = STARSH_SPATIAL_NON_GAUSSIAN_SIMD;
            //double noise = 0;
            STARSH_ssdata *ssdata;
            //int starsh_ssdata_init(STARSH_ssdata **data, STARSH_int count, int ndim,
            /*double *point, double beta, double nu, double noise, double sigma)*/
            info = starsh_ssdata_init(&ssdata, _M, ndim, point, initial_theta[0],
                                      initial_theta[1], noise, 1);
            /*info = starsh_application((void **)&ssdata, &kernel, _M, dtype,*/
            /*STARSH_SPATIAL, kernel_type, STARSH_SPATIAL_NDIM, ndim,*/
            /*STARSH_SPATIAL_BETA, initial_theta[1], STARSH_SPATIAL_NU, initial_theta[2],*/
            /*STARSH_SPATIAL_NOISE, noise,*/
            /*STARSH_SPATIAL_PLACE, STARSH_PARTICLES_OBSOLETE1, //1,*/
            /*STARSH_SPATIAL_SIGMA, initial_theta[0],*/
            /*0);*/
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_init\n");
            }
            /*int starsh_ssdata_get_kernel(STARSH_kernel **kernel, STARSH_ssdata *data,*/
            /*enum STARSH_SPATIAL_KERNEL type)*/
            info = starsh_ssdata_get_kernel(&kernel, ssdata, kernel_type);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_get_kernel\n");
            }


            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            //    printf("STARSH problem was succesfully generated\n");
            //    starsh_problem_info(problem);
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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

            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS_BIVARIATE_POINT) {
            //   printf("Geostat");
            //double initial_theta[3] = {1, 0.1, 0.5};

            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            //srand(0);  //no need in case og real dataset
            int kernel_type = STARSH_SPATIAL_PARSIMONIOUS_SIMD;
            //double noise = 0;
            STARSH_ssdata *ssdata;
            info = starsh_ssdata_init_parsimonious(&ssdata, _M, ndim, point, initial_theta[0],
                                                   initial_theta[1], initial_theta[2], initial_theta[3],
                                                   initial_theta[4], initial_theta[5], noise);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_init\n");
            }
            info = starsh_ssdata_get_kernel(&kernel, ssdata, kernel_type);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_get_kernel\n");
            }
            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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
            *mpiF = F;
        } else if (probtype == HICMA_STARSH_PROB_GEOSTAT_PARSIMONIOUS2_BIVARIATE_POINT) {
            //   printf("Geostat");
            //double initial_theta[3] = {1, 0.1, 0.5};

            char symm = 'S', dtype = 'd';
            STARSH_int shape[2] = {_M, _M};
            int info;
            //srand(0);  //no need in case og real dataset
            int kernel_type = STARSH_SPATIAL_PARSIMONIOUS2_SIMD;
            //double noise = 0;
            STARSH_ssdata *ssdata;
            info = starsh_ssdata_init_parsimonious(&ssdata, _M, ndim, point, initial_theta[0],
                                                   initial_theta[1], initial_theta[2], initial_theta[3],
                                                   initial_theta[4], initial_theta[5], noise);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_init\n");
            }
            info = starsh_ssdata_get_kernel(&kernel, ssdata, kernel_type);
            if (info != 0) {
                printf("wrong parameters for starsh_ssdata_get_kernel\n");
            }

            STARSH_problem *problem;
            info = starsh_problem_new(&problem, ndim, shape, symm, dtype, ssdata, ssdata, kernel,
                                      "Spatial Statistics example");
            if (info != 0) {
                printf("Error in starsh problem\n");
                exit(info);
            }
            STARSH_cluster *cluster;
            info = starsh_cluster_new_plain(&cluster, ssdata, _M, block_size);
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
            *mpiF = F;
        } else {
            fprintf(stderr, "Unknown type of STARS-H problem:%d. Exiting...\n", probtype);
        }
        //starsh_problem_info(P);
        //starsh_cluster_info(C);
    }
}

void HICMA_zgenerate_problem(
        int probtype,
        char sym,
        double decay,
        int _M,
        int _nb,
        int _mt,
        int _nt,
        HICMA_problem_t *hicma_problem
) {
    int diag_dense = 1;
    int initial_maxrank, final_maxrank;
    double initial_avgrank, final_avgrank;
#define ndescs  3
    double *lrarrays[ndescs] = {NULL, NULL, NULL};
    int num_descs = 1;
    HICMA_desc_t *descs[ndescs] = {NULL, NULL, NULL};
    HICMA_desc_t *descsU[ndescs] = {NULL, NULL, NULL};
    HICMA_desc_t *descsV[ndescs] = {NULL, NULL, NULL};
    HICMA_desc_t *descsD[ndescs] = {NULL, NULL, NULL};
    HICMA_desc_t *descsrk[ndescs] = {NULL, NULL, NULL};

    __generate(probtype, sym, decay, _M, _nb, _mt, _nt, num_descs, descs, descsU, descsV, descsD, descsrk, lrarrays,
               &initial_maxrank, &initial_avgrank, diag_dense, &(hicma_problem->starsh_format), hicma_problem->theta,
               hicma_problem->ndim, hicma_problem->beta, hicma_problem->nu, hicma_problem->noise,
               hicma_problem->kernel_type, hicma_problem->point, hicma_problem);
    HICMA_set_starsh_format(hicma_problem->starsh_format);

}
