/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file hicma_kernels.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon elementary kernels enum
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 */
#ifndef _HICMA_CHAM_KERNELS_H_
#define _HICMA_CHAM_KERNELS_H_

/**
 * Used to apply operations on specific kernels
 */
typedef enum hicma_kernel_e {

    HICMA_GEMM,

} HICMA_kernel_t;

#endif
