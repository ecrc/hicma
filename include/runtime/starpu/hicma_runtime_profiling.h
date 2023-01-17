/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_profiling.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU profiling and kernel locality header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _PROFILING_H_
#define _PROFILING_H_

#ifdef CHAMELEON_ENABLE_PRUNING_STATS
extern unsigned long RUNTIME_total_tasks;
extern unsigned long RUNTIME_exec_tasks;
extern unsigned long RUNTIME_comm_tasks;
extern unsigned long RUNTIME_changed_tasks;
#endif

typedef struct measure_s {
    double sum;
    double sum2;
    long   n;
} measure_t;

#endif
