# Release 1.0.0

## Features

* Removed Chameleon dependency completely.

# Release 0.1.4

## Features

* LU factorization (hicma_zgetrf) on double complex matrices stored in tile low-rank (TLR) format. ![Link to branch.](https://github.com/ecrc/hicma/tree/zgetrf)

# Release 0.1.3

## Features
* Support for a new application that simulates the 3D unstructured mesh deformation of a population of the novel coronaviruses (i.e.,  SARS-CoV-2). This relies on the virus geometry extracted from the Protein Data Bank (PDB) codenamed PDBID 6VXX available at (https://www.rcsb.org/structure/6VXX). The corresponding basis function and data geometry are available at STARS-H (https://github.com/ecrc/stars-h).

# Release 0.1.2

## Features
* New matrix kernel (3D exponential) for statistics application from STARS-H (https://github.com/ecrc/stars-h) is supported.
* Sequential operations on low-rank tiles are moved into a standalone library called HCORE (https://github.com/ecrc/hcore).
* Number of flops is reported. See the driver routines under timing folder. 
* Multiplication order of tiles in each inner product performed by HiCMA_GEMM can be changed via adding --reorderinnerproducts flag while calling timing/time_zgemm_tile. This flag enables sorting tiles according to increasing sum of ranks of tile pairs that which will be multiplied.
* More tests in timing/CMakeLists.txt 

## Bug Fixes and Other Changes
* Set minimum number of singular vectors in HCORE_GEMM to 1 instead of 2 

# Release 0.1.1

## Features
* Triangular solve (left) for low rank A and B/X
* Triangular solve (left, lower) for low rank A and full rank B/X  
* Testing routine for POSV and TRSM
* HCORE routine for C=C+AB where A is low rank and B and C are full rank 
* Separate routines for dense and low rank matrix generation 

# Release 0.1.0

## Features
* Matrix-Matrix Multiplication
* Cholesky Factorization
* Double Precision
* Task-based Programming Models
* Shared and Distributed-Memory Environments
* Support for StarPU Dynamic Runtime Systems
* Testing Suite and Examples

