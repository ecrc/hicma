/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

#include <hicma_init.h>
#include <hicma_common.h>

int HICMA_set_print_index(){
    hicma_context.print_index = 1;
    return 0;
}
int HICMA_unset_print_index(){
    hicma_context.print_index = 0;
    return 0;
}
int HICMA_set_print_index_end(){
    hicma_context.print_index_end = 1;
    return 0;
}
int HICMA_unset_print_index_end(){
    hicma_context.print_index_end = 0;
    return 0;
}
int HICMA_set_use_fast_hcore_gemm(){
    hicma_context.use_fast_hcore_gemm = 1;
    return 0;
}
int HICMA_get_use_fast_hcore_gemm(){
    return hicma_context.use_fast_hcore_gemm;
}
int HICMA_unset_use_fast_hcore_gemm(){
    hicma_context.use_fast_hcore_gemm = 0;
    return 0;
}
int HICMA_get_always_fixed_rank(){
    return hicma_context.global_always_fixed_rank;
}
int HICMA_get_fixed_rank(){
    return hicma_context.global_fixed_rank;
}
int HICMA_set_fixed_rank(int rank){
    hicma_context.global_fixed_rank = rank;
    return 0;
}
int HICMA_get_print_index(){
    return hicma_context.print_index;
}
int HICMA_get_print_index_end(){
    return hicma_context.print_index_end;
}
int HICMA_get_print_mat(){
    return hicma_context.print_mat;
}
int HICMA_set_print_mat(){
    hicma_context.print_mat = 1;
    return 0;
}
int HICMA_set_starsh_format(STARSH_blrf *starsh_format){
    hicma_context.starsh_format = starsh_format; 
    return 0;
}
STARSH_blrf * HICMA_get_starsh_format(){
    return hicma_context.starsh_format; 
}


/**
 *
 * @ingroup Control
 *
 *  HICMA_Init - Initialize HICMA.
 *
 ******************************************************************************
 *
 * @param[in] cores
 *          Number of cores to use.
 *
 * @param[in] gpus
 *          Number of cuda devices to use.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Init(int cores, int gpus)
{
    hicma_context = hicma_context_default;
    return HICMA_InitPar(cores, gpus, -1);
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_InitPar - Initialize HICMA.
 *
 ******************************************************************************
 *
 * @param[in] ncpus
 *          Number of cores to use.
 *
 * @param[in] ncudas
 *          Number of cuda devices to use.
 *
 * @param[in] nthreads_per_worker
 *          Number of threads per worker (cpu, cuda device).
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_InitPar(int ncpus, int ncudas, int nthreads_per_worker)
{
    HICMA_context_t *hicma;

    /* Create context and insert in the context map */
    hicma = hicma_context_create();
    if (hicma == NULL) {
        hicma_fatal_error("HiCMA_Init", "hicma_context_create() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }

#if defined(HICMA_USE_MPI)
    #  if defined(CHAMELEON_SIMULATION)
    /* Assuming that we don't initialize MPI ourself (which SMPI doesn't support anyway) */
    hicma->mpi_outer_init = 1;
#  else
    {
      int flag = 0, provided = 0;
      MPI_Initialized( &flag );
      hicma->mpi_outer_init = flag;
      if ( !flag ) {
          MPI_Init_thread( NULL, NULL, MPI_THREAD_MULTIPLE, &provided );
      }
    }
#  endif
#endif

    HICMA_RUNTIME_init( hicma, ncpus, ncudas, nthreads_per_worker );

#if defined(HICMA_USE_MPI)
    hicma->my_mpi_rank   = HICMA_RUNTIME_comm_rank( hicma );
    hicma->mpi_comm_size = HICMA_RUNTIME_comm_size( hicma );
#endif

    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Finalize - Finalize HICMA.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Finalize(void)
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Finalize()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    HICMA_RUNTIME_flush();
#  if !defined(CHAMELEON_SIMULATION)
    HICMA_RUNTIME_barrier(hicma);
#  endif
    HICMA_RUNTIME_finalize( hicma );

#if defined(HICMA_USE_MPI)
    if (!hicma->mpi_outer_init)
        MPI_Finalize();
#endif

    hicma_context_destroy();
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Pause - Suspend HICMA runtime to poll for new tasks.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Pause(void)
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Pause()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    HICMA_RUNTIME_pause(hicma);
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Resume - Symmetrical call to HICMA_Pause,
 *  used to resume the workers polling for new tasks.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Resume(void)
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Resume()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    HICMA_RUNTIME_resume(hicma);
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Distributed_start - Prepare the distributed processes for computation
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Distributed_start(void)
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Finalize()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    HICMA_RUNTIME_barrier (hicma);
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Distributed_stop - Clean the distributed processes after computation
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Distributed_stop(void)
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Finalize()", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    HICMA_RUNTIME_barrier (hicma);
    return HICMA_SUCCESS;
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Comm_size - Return the size of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The size of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int HICMA_Comm_size()
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Comm_size()", "HiCMA not initialized");
        return -1;
    }

    return HICMA_RUNTIME_comm_size( hicma );
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_Comm_rank - Return the rank of the distributed computation
 *
 ******************************************************************************
 *
 * @retval The rank of the distributed computation
 * @retval -1 if context not initialized
 *
 */
int HICMA_Comm_rank()
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Comm_rank()", "HiCMA not initialized");
        return -1;
    }

    return HICMA_RUNTIME_comm_rank( hicma );
}

/**
 *
 * @ingroup Control
 *
 *  HICMA_GetThreadNbr - Return the number of CPU workers initialized by the
 *  runtime
 *
 ******************************************************************************
 *
 * @return
 *          \retval The number of CPU workers started
 *
 */
int HICMA_GetThreadNbr( )
{
    HICMA_context_t *hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_GetThreadNbr()", "HiCMA not initialized");
        return -1;
    }

    return HICMA_RUNTIME_thread_size( hicma );
}
