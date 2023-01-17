/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_workspace.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU workspace header
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @date 2011-06-01
 *
 */
#ifndef _HICMA_CHAM_STARPU_WORKSPACE_H_
#define _HICMA_CHAM_STARPU_WORKSPACE_H_

/*
 * Allocate workspace in host memory: CPU for any worker
 * or allocate workspace in worker's memory: main memory for cpu workers,
 * and embedded memory for CUDA devices.
 */
#define HICMA_HOST_MEM    0
#define HICMA_WORKER_MEM  1

struct hicma_starpu_ws_s {
    size_t size;
    int    memory_location;
    int    which_workers;
    void  *workspaces[STARPU_NMAXWORKERS];
};

typedef struct hicma_starpu_ws_s HICMA_starpu_ws_t;

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. HICMA_CUDA|HICMA_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in worker's memory (HICMA_HOST_MEM or HICMA_WORKER_MEM). This function
 * returns 0 upon successful completion.
 */
int   HICMA_RUNTIME_starpu_ws_alloc   ( HICMA_starpu_ws_t **workspace, size_t size, int which_workers, int memory_location);
int   HICMA_RUNTIME_starpu_ws_free    ( HICMA_starpu_ws_t  *workspace);
void *HICMA_RUNTIME_starpu_ws_getlocal( HICMA_starpu_ws_t  *workspace);

#endif /* _HICMA_CHAM_STARPU_WORKSPACE_H_ */
