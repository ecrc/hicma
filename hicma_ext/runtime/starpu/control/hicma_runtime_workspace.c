/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_workspace.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU workspaces routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @date 2011-06-01
 *
 */
#include "runtime/starpu/hicma_starpu.h"

static void HICMA_RUNTIME_allocate_workspace_on_workers(void *arg)
{
    struct hicma_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;

    int id = starpu_worker_get_id();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    type = starpu_worker_get_type(id);
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;

        if (memory_location == HICMA_HOST_MEM)
        {
            /* Use pinned memory because the kernel is very likely
             * to transfer these data between the CPU and the GPU.
             * */
            cudaMallocHost(&workspace->workspaces[id], workspace->size);
        }
        else {
            /* Allocate on the device */
            cudaMalloc(&workspace->workspaces[id], workspace->size);
        }
    }
    else
#endif
    {
        /* This buffer should only be used within the CPU kernel, so
         * there is no point in using pinned memory here. */
        workspace->workspaces[id] = malloc(workspace->size);
    }

    assert(workspace->workspaces[id]);
}


static void HICMA_RUNTIME_free_workspace_on_workers(void *arg)
{
    struct hicma_starpu_ws_s *workspace = arg;
    enum starpu_worker_archtype type = 0;
    (void)type;
    int id = starpu_worker_get_id();

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    type = starpu_worker_get_type(id);
    if (type == STARPU_CUDA_WORKER)
    {
        int memory_location = workspace->memory_location;

        if (memory_location == HICMA_HOST_MEM)
        {
            cudaFreeHost(workspace->workspaces[id]);
        }
        else {
            cudaFree(workspace->workspaces[id]);
        }
    }
    else
#endif
    {
        free(workspace->workspaces[id]);
    }

    workspace->workspaces[id] = NULL;
}

/*
 * This function creates a workspace on each type of worker in "which_workers"
 * (eg. HICMA_CUDA|HICMA_CPU for all CPU and GPU workers).  The
 * memory_location argument indicates whether this should be a buffer in host
 * memory or in GPU memory (HICMA_HOST_MEM or HICMA_GPU_MEM). This function
 * returns 0 upon successful completion.:
 */
int HICMA_RUNTIME_starpu_ws_alloc(HICMA_starpu_ws_t **workspace,
                            size_t size, int which_workers, int memory_location)
{
    if (!workspace)
        return -EINVAL;

    struct hicma_starpu_ws_s *descr = calloc(1, sizeof(struct hicma_starpu_ws_s));

    *workspace = descr;

    if (!descr)
        return -ENOMEM;

    descr->size = size;
    descr->memory_location = memory_location;
    descr->which_workers = which_workers;

    starpu_execute_on_each_worker(HICMA_RUNTIME_allocate_workspace_on_workers, descr, which_workers);

    return 0;
}

int HICMA_RUNTIME_starpu_ws_free(HICMA_starpu_ws_t *workspace)
{
    if (!workspace)
        return -EINVAL;

    starpu_execute_on_each_worker(HICMA_RUNTIME_free_workspace_on_workers, workspace, workspace->which_workers);

    free(workspace);

    return 0;
}

void *HICMA_RUNTIME_starpu_ws_getlocal(HICMA_starpu_ws_t *workspace)
{
    struct hicma_starpu_ws_s *descr = workspace;
    int id = starpu_worker_get_id();
    return descr->workspaces[id];
}
