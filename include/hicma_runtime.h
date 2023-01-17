/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file hicma_runtime.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2017 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon runtimes API
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @date 2011-06-01
 *
 */
#ifndef _HICMA_CHAM_RUNTIME_H_
#define _HICMA_CHAM_RUNTIME_H_

#include "hicma_config.h"
#include <hicma_struct.h>

BEGIN_C_DECLS

/**
 * @name RUNTIME Context functions
 * @{
 *    These functions manipulate the context data structure before the runtime is
 *    started and after the runtime is stopped.
 */

/**
 * @brief Create the runtime specific options in the context before starting it
 *
 * @param[in,out] ctxt
 *            The Chameleon context to initialize for the runtime.
 */
void
HICMA_RUNTIME_context_create( HICMA_context_t *ctxt );

/**
 * @brief Destroy the specific options in the context after this last one has
 * been stop.
 *
 * @param[in,out] ctxt
 *            The Chameleon context in which the runtime specific options must
 *            be destroyed
 */
void
HICMA_RUNTIME_context_destroy( HICMA_context_t *ctxt );

/**
 * @brief Enable a global option of the runtime.
 * @warning Should be called only by HICMA_Enable()
 *
 * @param[in] option
 *            @arg HICMA_PROFILING_MODE: start the profiling mode of the runtime.
 */
void
HICMA_RUNTIME_enable( HICMA_enum option );

/**
 * @brief Disable a global option of the runtime.
 * @warning Should be called only by HICMA_Disable()
 *
 * @param[in] option
 *            @arg HICMA_PROFILING_MODE: stop the profiling mode of the runtime.
 */
void
HICMA_RUNTIME_disable( HICMA_enum option );

/**
 * @}
 *
 * @name RUNTIME Control functions
 * @{
 *   These functions control the global behavior of the runtime.
 */

/**
 * @brief Initialize the scheduler with the given parameters
 *
 * @param[in,out] ctxt
 *            The Chameleon context in which to initialize the runtime support.
 *
 * @param[in] ncpus
 *            Defines the total number of cores given to the runtime per
 *            node. Including cuda and communication workers for runtimes that
 *            dedicates cores to this. ncpus > 0, or -1 to target the whole machine.
 *
 * @param[in] ncudas
 *            Defines the number of CUDA devices used by node. If the runtime
 *            uses one core dedicated to each CUDA device, they will be taken
 *            from ncpus. If ncpus > 0, ncudas < ncpus. -1 to target all the
 *            CUDA devices available.
 *
 * @param[in] nthreads_per_worker
 *            Defines the number of threads per worker when multi-threaded tasks
 *            are enabled.  This is used to exploit parallel BLAS kernels, and
 *            defines a better binding of the workers.
 *            -1 to disable, or > 0 to enable.
 *
 * @retval -1 on failure to initialize the runtime.
 * @retval >0 on success to initialize the runtime.
 *
 */
int
HICMA_RUNTIME_init( HICMA_context_t *ctxt,
              int ncpus,
              int ncudas,
              int nthreads_per_worker );

/**
 * @brief Finalize the scheduler used for the computations.
 *
 * @param[in,out] ctxt
 *            The Chameleon context for which the runtime system must be shut down.
 */
void
HICMA_RUNTIME_finalize( HICMA_context_t *ctxt );

/**
 * @brief Suspend the processing of new tasks submitted to the runtime system.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the suspension must be made.
 */
void
HICMA_RUNTIME_pause( HICMA_context_t *ctxt );

/**
 * @brief Resume the processing of new tasks submitted to the runtime system.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the execution must be resumed.
 */
void
HICMA_RUNTIME_resume( HICMA_context_t *ctxt );

/**
 * @brief Wait for completion of all tasks submitted to the runtime.
 *
 * @param[in] ctxt
 *            The Chameleon context in which the task completion is performed.
 */
void
HICMA_RUNTIME_barrier( HICMA_context_t *ctxt );

/**
 * @brief Show the progress of the computations when enabled.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the context needs to be printed.
 */
void
HICMA_RUNTIME_progress( HICMA_context_t *ctxt );

/**
 * @brief Get the rank of the current worker for the runtime.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the thread rank is asked.
 *
 * @retval The rank of the current thread in the runtime.
 */
int
HICMA_RUNTIME_thread_rank( HICMA_context_t *ctxt );

/**
 * @brief Get the number of CPU workers of the runtime.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the number of workers is requested
 *
 * @retval The number of threads currently used by the runtime.
 */
int
HICMA_RUNTIME_thread_size( HICMA_context_t *ctxt );

/**
 * @brief Get the MPI comm rank of the current process related to the runtime.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the rank is asked.
 *
 * @retval The rank of the process in the communicator known by the runtime.
 */
int
HICMA_RUNTIME_comm_rank( HICMA_context_t *ctxt );

/**
 * @brief Get the MPI comm size related to the runtime.
 *
 * @param[in] ctxt
 *            The Chameleon context for which the communicator size is asked.
 *
 * @retval The size of the communicator known by the runtime.
 */
int
HICMA_RUNTIME_comm_size( HICMA_context_t *ctxt );

/**
 * @brief Set the data tags bounds based on runtime tags limit.
 *
 * @param[in] user_tag_width
 *            Total number of bits used to defined the data tags (31 by default).
 *            Must be larger than 20.
 *
 * @param[in] user_tag_sep
 *            The number of bits dedicated to identify the pieces of data per
 *            descriptor.  (user_tag_width - user_tag_sep) defines the number of
 *            bits used to enumerate the descriptors.
 */
void
HICMA_RUNTIME_comm_set_tag_sizes( int user_tag_width,
                            int user_tag_sep );

/**
 * @}
 *
 * @name RUNTIME Asynchonous functionalities
 * @{
 *    These functions manage the sequences of tasks. A sequence is a subset of
 *    related tasks belonging to a same problem.
 */

/**
 * @brief Create a sequence structure associated to a given context.
 *
 * @param[in] ctxt
 *            The Chameleon context in which the sequence is created.
 *
 * @param[in,out] sequence
 *            On entry the allocated Chameleon sequence structure, and on exit
 *            the scheduler specifics for the sequence have been initialized.
 *
 * @retval HICMA_SUCCESS on success.
 * @retval HICMA_ERR_OUT_OF_RESOURCES, if the sequence could not be created.
 */
int
HICMA_RUNTIME_sequence_create( HICMA_context_t  *ctxt,
                         HICMA_sequence_t *sequence );

/**
 * @brief Destroy the sequence structure.
 *
 * @param[in] ctxt
 *            The Chameleon context to which the sequence belongs.
 *
 * @param[in,out] sequence
 *            On entry the sequence structure.
 *            On exit, the scheduler specifics of the sequence have been
 *            destroyed.
 *
 * @retval HICMA_SUCCESS on success.
 */
int
HICMA_RUNTIME_sequence_destroy( HICMA_context_t  *ctxt,
                          HICMA_sequence_t *sequence);

/**
 * @brief Wait for completion of all tasks in the given sequence.
 *
 * @param[in] ctxt
 *            The Chameleon context to which the sequence belongs.
 *
 * @param[in] sequence
 *            The sequence that gathers a set of tasks to complete.
 *
 * @retval HICMA_SUCCESS on success.
 * @retval The algorithm return code on failure.
 */
int
HICMA_RUNTIME_sequence_wait( HICMA_context_t  *ctxt,
                       HICMA_sequence_t *sequence );

/**
 * @brief Cancel all remaining tasks in the given sequence.
 *
 * @param[in] ctxt
 *            The Chameleon context to which the sequence belongs.
 *
 * @param[in,out] sequence
 *            The sequence that gathers the set of tasks to cancel.
 *            On exit, the sequence holds the error code of the algorithm and
 *            the request that failed in the sequence.
 *
 * @param[in,out] request
 *            The request that failed in the sequence.
 *            On exit, the return status is set to status.
 *
 * @param[in] status
 *            The return status of the task that failed in the request given
 *            above.
 */
void
HICMA_RUNTIME_sequence_flush( HICMA_context_t  *ctxt,
                        HICMA_sequence_t *sequence,
                        HICMA_request_t  *request,
                        int               status );

/**
 * @}
 *
 * @name RUNTIME Memory management
 * @{
 ***
 * @brief Allocate size bytes through the runtime memory management system if any, or with malloc otherwise.
 *
 * This function allows to allocate pinned memory if needed, or eventually not
 * perform the allocation in simulation mode.

 * @param[in] size
 *            The number of bytes to allocate.
 *
 * @return The pointer to allocated area of size bytes on success, NULL otherwise.
 */
void *
HICMA_RUNTIME_malloc( size_t size );

/**
 * @brief Free allocated memory through HICMA_RUNTIME_malloc() function call
 *
 * @param[in,out] ptr
 *            The ptr to free
 *
 * @param[in] size
 *            The size in bytes of the allocated area associated to ptr.
 */
void
HICMA_RUNTIME_free( void *ptr, size_t size );

/**
 * @}
 *
 * @name RUNTIME Descriptor functions
 * @{
 ***
 * @brief Initialize runtime specific data structure to a given descriptor.
 *
 * @param[in,out] desc
 *            The descriptor to initialize.
 */
void
HICMA_RUNTIME_desc_create( HICMA_desc_t *desc );

/**
 * @brief Finalize runtime specific data structure of a given descriptor.
 *
 * @param[in,out] desc
 *            The descriptor to finalize.
 */
void
HICMA_RUNTIME_desc_destroy( HICMA_desc_t *desc );

/**
 * @brief Acquire in main memory an up-to-date copy of the data described by the
 * descriptor for read-write access.
 *
 * The application must call this function prior to accessing registered data
 * from main memory outside tasks. HICMA_RUNTIME_desc_release() must be called once
 * the application does not need to access the data anymore.
 *
 * @param[in] desc
 *            The descriptor to acquire.
 *
 * @retval HICMA_SUCCESS on success
 */
int
HICMA_RUNTIME_desc_acquire( const HICMA_desc_t *desc );

/**
 * @brief Release the data described by the descriptor to be used by the runtime
 * tasks again.
 *
 * This function releases the data acquired by the application either by
 * HICMA_RUNTIME_desc_acquire() or by HICMA_RUNTIME_desc_acquire_async() to the runtime.
 *
 * @param[in] desc
 *            The descriptor to release.
 *
 * @retval HICMA_SUCCESS on success
 */
int
HICMA_RUNTIME_desc_release( const HICMA_desc_t *desc );

/**
 * @brief Flush all pieces of data from a descriptor.
 *
 * This function marks all pieces of data from the descriptor as unused in the future to:
 *   - cleanup the data from the distributed cache if present.
 *   - mark for eviction from the GPU if space is needed
 *   - move the data back to the main memory
 *
 * @param[in] desc
 *            The descriptor to flush.
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 */
void
HICMA_RUNTIME_desc_flush( const HICMA_desc_t     *desc,
                    const HICMA_sequence_t *sequence );

/**
 * @brief Flush all data submitted to the runtime systems from the distributed
 * cache, and/or mark for eviction from the GPU memory.
 *
 * This function flushes all data from the distributed cache of the runtime system.
 */
void
HICMA_RUNTIME_flush( );

/**
 * @brief Flush a single piece of data.
 *
 * This function marks a piece of data as unused in the future to:
 *   - cleanup the data from the distributed cache if present.
 *   - mark for eviction from the GPU if space is needed
 *   - move the data back to the main memory
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 *
 * @param[in] A
 *            The descriptor to which the piece of data belongs.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 */
void
HICMA_RUNTIME_data_flush( const HICMA_sequence_t *sequence,
                    const HICMA_desc_t *A, int Am, int An );

/**
 * @brief Migrate a single piece of data.
 *
 * This function migrate a piece of data from its original rank to the new_rank
 * and changes its ownership.
 *
 * @param[in] sequence
 *            The sequence in which the data is used.
 *
 * @param[in] A
 *            The descriptor to which the piece of data belongs.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 *
 * @param[in] new_rank
 *            The new_rank on which to migrate the data
 */
#if defined(CHAMELEON_USE_MIGRATE)
void
HICMA_RUNTIME_data_migrate( const HICMA_sequence_t *sequence,
                      const HICMA_desc_t *A, int Am, int An, int new_rank );
#else
static inline void
HICMA_RUNTIME_data_migrate( const HICMA_sequence_t *sequence,
                      const HICMA_desc_t *A, int Am, int An, int new_rank ) {
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
}
#endif

/**
 * @brief Get the pointer to the data or the runtime handler associated to the
 * piece of data (m, n) in desc.
 *
 * @param[in] A
 *            The descriptor to which belongs the piece of data.
 *
 * @param[in] Am
 *            The row coordinate of the piece of data in the matrix
 *
 * @param[in] An
 *            The column coordinate of the piece of data in the matrix
 *
 * @retval The runtime handler address of the piece of data.
 */
void *
HICMA_RUNTIME_data_getaddr( const HICMA_desc_t *A, int Am, int An );

/**
 * @}
 *
 * @name RUNTIME Insert task options management
 * @{
 */
void  HICMA_RUNTIME_options_init     (HICMA_option_t*, HICMA_context_t*, HICMA_sequence_t*, HICMA_request_t*);
void  HICMA_RUNTIME_options_finalize (HICMA_option_t*, HICMA_context_t *);
int   HICMA_RUNTIME_options_ws_alloc (HICMA_option_t*, size_t, size_t);
int   HICMA_RUNTIME_options_ws_free  (HICMA_option_t*);
/* int   RUNTIME_options_ws_gethost   (HICMA_option_t*); */
/* int   RUNTIME_options_ws_getdevice (HICMA_option_t*); */

/**
 * @}
 *
 * @name RUNTIME Kernel locality management
 * @{
 */
void RUNTIME_zlocality_allrestore ();
void RUNTIME_clocality_allrestore ();
void RUNTIME_dlocality_allrestore ();
void RUNTIME_slocality_allrestore ();

void RUNTIME_zlocality_allrestrict(uint32_t);
void HICMA_RUNTIME_zlocality_onerestrict(HICMA_kernel_t, uint32_t);
void HICMA_RUNTIME_zlocality_onerestore (HICMA_kernel_t);

void RUNTIME_clocality_allrestrict(uint32_t);
void HICMA_RUNTIME_clocality_onerestrict(HICMA_kernel_t, uint32_t);
void HICMA_RUNTIME_clocality_onerestore (HICMA_kernel_t);

void RUNTIME_dlocality_allrestrict(uint32_t);
void HICMA_RUNTIME_dlocality_onerestrict(HICMA_kernel_t, uint32_t);
void HICMA_RUNTIME_dlocality_onerestore (HICMA_kernel_t);

void RUNTIME_slocality_allrestrict(uint32_t);
void HICMA_RUNTIME_slocality_onerestrict(HICMA_kernel_t, uint32_t);
void HICMA_RUNTIME_slocality_onerestore (HICMA_kernel_t);

/**
 * @}
 *
 * @name RUNTIME Profiling
 * @{
 */
void   HICMA_RUNTIME_schedprofile_display ();
void   HICMA_RUNTIME_kernelprofile_display();
double HICMA_RUNTIME_get_time();

void  HICMA_RUNTIME_iteration_push     (HICMA_context_t*, unsigned long iteration);
void  HICMA_RUNTIME_iteration_pop      (HICMA_context_t*);

void HICMA_RUNTIME_start_profiling();
void HICMA_RUNTIME_stop_profiling();

void HICMA_RUNTIME_start_stats();
void HICMA_RUNTIME_stop_stats();

void RUNTIME_zdisplay_allprofile ();
void HICMA_RUNTIME_zdisplay_oneprofile (HICMA_kernel_t);
void RUNTIME_cdisplay_allprofile ();
void HICMA_RUNTIME_cdisplay_oneprofile (HICMA_kernel_t);
void RUNTIME_ddisplay_allprofile ();
void HICMA_RUNTIME_ddisplay_oneprofile (HICMA_kernel_t);
void RUNTIME_sdisplay_allprofile ();
void HICMA_RUNTIME_sdisplay_oneprofile (HICMA_kernel_t);

/**
 * @}
 */

END_C_DECLS

#endif /* _HICMA_RUNTIME_H_ */
