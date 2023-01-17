/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

#include <stdio.h>
#include <starsh_core.h>
// Add this header, wait it to be fixed by Aleks
//#include <starsh_core/core/check.h>

// ~/stars-h-core-dev/examples/c/applications/electrodynamics.c
int main(int argc, char **argv)
{
    // Space dimensionality for particles
    starsh_int space_ndim = 3;
    // Number of particles
    starsh_int particles_count = 10000;
    // Spatial coordinates in double precision
    starsh_dparticles_t *particles;
    // Generate coordinates in a uniform way
    starsh_dparticles_uniform(&particles, space_ndim,
            particles_count);
    // Sort coordinates with Z-order inplace
    starsh_dparticles_zsort_inplace(particles, space_ndim);
    // Electrodynamics object
    starsh_delectrodynamics_t *electrodynamics;
    // Parameters for geospatial statistics
    double wave_k = 10.0; // Wave number
    // Use coordinates to init electrodynamics object, which holds a copy
    // of provided coordinates
    starsh_delectrodynamics_constructor2(&electrodynamics, space_ndim,
            particles, wave_k);
    // Free coordinates
    starsh_dparticles_destructor(&particles, space_ndim);
    // Select complex matrix kernel for 3-dimensional space
    // Cast type explicitly to silence incompatible types warning
    starsh_zkernel_t *kernel =
        (starsh_zkernel_t *)starsh_delectrodynamics_zkernel_exp_3d;
    // Problem object, that corresponds to a matrix in a matrix-free form
    starsh_zproblem_t *problem;
    // Shape of corresponding matrix (ndim is a dimensionality of matrix)
    starsh_int problem_ndim = 2;
    starsh_int problem_shape[2] = {particles_count, particles_count};
    // Create starsh problem, a matrix in a matrix-free form
    // `spatial` is NOT copied, but stored as a reference/pointer
    starsh_zproblem_constructor(&problem, problem_ndim, problem_shape,
            electrodynamics, electrodynamics, kernel);
    // Generate a submatrix on first 10 rows and first 10 columns
    starsh_int nrows = 10;
    starsh_int ncols = nrows;
    starsh_int *index = malloc(nrows * sizeof(*index));
    for(starsh_int i = 0; i < nrows; ++i)
        index[i] = i;
    double _Complex *submatrix = malloc(nrows * ncols * sizeof(*submatrix));
    starsh_int ld = nrows; // Leading dimension of submatrix
    char trans = 'N'; // no transposition
    // Generate submatrix
    starsh_zproblem_call_kernel(problem, trans, nrows, index, nrows, index,
            submatrix, nrows);
    // Now we generate low-rank off-diagonal submatrix in the left bottom corner
    starsh_int *index2 = malloc(nrows * ncols * sizeof(*index2));
    for(starsh_int i = 0; i < nrows; ++i)
        index2[i] = particles_count - nrows + i;
    starsh_zproblem_call_kernel(problem, trans, nrows, index2, ncols, index,
            submatrix, nrows);
    // Approximate it now
    starsh_int maxrank = 10;
    starsh_int oversample = 5;
    double _Complex *U = malloc(nrows * maxrank * sizeof(*U));
    double _Complex *V = malloc(ncols * maxrank * sizeof(*V));
    starsh_int rank = -1;
    starsh_int lwork = -1;
    lapack_int iseed[4] = {0, 0, 0, 1};
    double tol = 1e-7;
    char norm = 'F';
    // Compute optimal workspace size
    starsh_zcore_rsdd_work(nrows, ncols, NULL, nrows, NULL, nrows, NULL, ncols,
            NULL, maxrank, oversample, NULL, norm, tol, NULL, &lwork);
    // Now allocate workspace and run approximate
    double _Complex *work = malloc(lwork * sizeof(*work));
    starsh_zcore_rsdd_work(nrows, ncols, submatrix, nrows, U, nrows, V, ncols,
            &rank, maxrank, oversample, iseed, norm, tol, work, &lwork);
    // Check accuracy
    double err_norm = -1;
    starsh_zcore_check(nrows, ncols, submatrix, nrows, rank, U, nrows, V, ncols,
            &err_norm);
    printf("Accuracy: %e\nRank: %d\n", err_norm, rank);
    // Free problem
    starsh_zproblem_destructor(&problem);
    // Free electrodynamics object
    starsh_delectrodynamics_destructor(&electrodynamics, space_ndim);
}
