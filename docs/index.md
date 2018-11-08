# What is HiCMA?

The Hierarchical Computations on Manycore Architectures (HiCMA) library aims to redesign existing dense linear algebra libraries to exploit the data sparsity of the matrix operator. Data sparse matrices arise in many scientific problems (e.g., in statistics-based weather forecasting, seismic imaging, and materials science applications) and are characterized by low-rank off-diagonal tile structure. Numerical low-rank approximations have demonstrated attractive theoretical bounds,
both in memory footprint and arithmetic complexity. The core idea of HiCMA is to develop fast linear algebra computations operating on the underlying tile low-rank data format, while satisfying a specified numerical accuracy and leveraging performance from massively parallel hardware architectures.

# Features of HiCMA 0.1.1

* Matrix-Matrix Multiplication
* Cholesky Factorization (Solve will be available soon)
* Double Precision
* Task-based Programming Models
* Shared and Distributed-Memory Environments
* Support for StarPU Dynamic Runtime Systems
* Testing Suite and Examples

# Installation

Installation requires `CMake` of version 3.2.3 at least. To build HiCMA,
follow these instructions:

1.  Get HiCMA from git repository

        git clone git@github.com:ecrc/hicma


2.  Go into hicma folder

        cd hicma

3.  Get submodules using git as follows. The submodules Chameleon and STARS-H should be compiled and the `PKG_CONFIG_PATH` should be set.

        git submodule update --init --recursive

4.  Create build directory and go there

        mkdir build && cd build

5.  Use CMake to get all the dependencies

        cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install/  -DHICMA_USE_MPI=ON

6.  Build HiCMA

        make -j

7.  Build local documentation (optional)

        make docs

8.  Install HiCMA

        make install

9. Add line

        export PKG_CONFIG_PATH=/path/to/install:$PKG_CONFIG_PATH

    to your .bashrc file to use HiCMA as a library.

Now you can use `pkg-config` executable to collect compiler and linker flags for HiCMA or run the binaries under `/path/to/install/timing/`.

There are `scripts/build.sh` and `scripts/build-nompi.sh` scripts in the repository to build the whole software stack.

# Quick Start

## Cholesky Factorization 
    timing/time_zpotrf_tile.c performs the following operations:
    1. Creates a matrix in tile low-rank (TLR) format.
    2. Performs Cholesky factorization on the TLR matrix.
    3. Checks solution according to LAPACK_dpotrf() if checking is enabled.

## General Matrix Multiply
    timing/time_zgemm_tile.c performs the following operations:
    1. Creates three matrices in tile low-rank (TLR) format.
    2. Performs matrix multiplication.
    3. Checks solution according to cblas_gemm() if checking is enabled.
