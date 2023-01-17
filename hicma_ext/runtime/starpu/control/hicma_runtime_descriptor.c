/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file runtime_descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon StarPU descriptor routines
 *
 * @version 1.0.0
 * @author Cedric Augonnet
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 */
#include <stdlib.h>
#include <unistd.h>
#include "runtime/starpu/hicma_starpu.h"

/**
 *  Set the tag sizes
 */
#if defined(HICMA_USE_MPI)

/* Take 24 bits for the tile id, and 7 bits for descriptor id.
 These values can be changed through the call HICMA_user_tag_size(int tag_width, int tag_sep) */
#define TAG_WIDTH_MIN 20
static int tag_width = 31;
static int tag_sep   = 24;
static int _tag_mpi_initialized_ = 0;

static inline int
chameleon_starpu_tag_init( int user_tag_width,
                           int user_tag_sep )
{
    if (!_tag_mpi_initialized_) {
        int ok = 0;
        uintptr_t tag_ub;

        tag_width = user_tag_width;
        tag_sep   = user_tag_sep;

#if defined(HAVE_STARPU_MPI_COMM_GET_ATTR)
        int64_t *tag_ub_p = NULL;
        starpu_mpi_comm_get_attr(MPI_COMM_WORLD, STARPU_MPI_TAG_UB, &tag_ub_p, &ok);
        tag_ub = *tag_ub_p;
#else
        int *tag_ub_p = NULL;
        MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_TAG_UB, &tag_ub_p, &ok);
        tag_ub = *tag_ub_p;
#endif

        if ( !ok ) {
            hicma_error("HiCMA_cham_RUNTIME_desc_create", "MPI_TAG_UB not known by StarPU");
        }

        while ( ((uintptr_t)((1UL<<tag_width) - 1) > tag_ub ) &&
                (tag_width >= TAG_WIDTH_MIN) )
        {
            tag_width--;
            tag_sep--;
        }

        if ( tag_width < TAG_WIDTH_MIN ) {
            hicma_error("HiCMA_cham_RUNTIME_desc_create", "MPI_TAG_UB too small to identify all the data");
            return HICMA_ERR_OUT_OF_RESOURCES;
        }

        _tag_mpi_initialized_ = 1;
        return HICMA_SUCCESS;
    }
    else {
        return HICMA_ERR_REINITIALIZED;
    }
}


#ifndef HAVE_STARPU_MPI_DATA_REGISTER
#define starpu_mpi_data_register( handle_, tag_, owner_ )       \
    do {                                                        \
        starpu_data_set_rank( (handle_), (owner_) );            \
        starpu_data_set_tag( (handle_), (tag_) );               \
    } while(0)
#endif

#endif

void HICMA_RUNTIME_comm_set_tag_sizes(int user_tag_width,
                                           int user_tag_sep) {
#if defined(HICMA_USE_MPI)
    int rc;
    rc = chameleon_starpu_tag_init( user_tag_width, user_tag_sep );
    if ( rc != HICMA_SUCCESS ) {
        hicma_error("RUNTIME_user_tag_size",
                    "must be called before creating any Hicma descriptor with HICMA_Desc_create(). The tag sizes will not be modified.");
    }
#endif
    (void) user_tag_width;
    (void) user_tag_sep;
}

/**
 *  Malloc/Free of the data
 */
#ifdef STARPU_MALLOC_SIMULATION_FOLDED
#define FOLDED STARPU_MALLOC_SIMULATION_FOLDED
#else
#define FOLDED 0
#endif

void *HICMA_RUNTIME_malloc(size_t size) {
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(HICMA_USE_MPI)
    return (void*) 1;
#else
    void *ptr;

    if (starpu_malloc_flags(&ptr, size, STARPU_MALLOC_PINNED | FOLDED | STARPU_MALLOC_COUNT) != 0) {
        return NULL;
    }
    return ptr;
#endif
}

void HICMA_RUNTIME_free(void *ptr,
                             size_t size) {
#if defined(CHAMELEON_SIMULATION) && !defined(STARPU_MALLOC_SIMULATION_FOLDED) && !defined(HICMA_USE_MPI)
    (void)ptr; (void)size;
    return;
#else
    starpu_free_flags(ptr, size, STARPU_MALLOC_PINNED | FOLDED | STARPU_MALLOC_COUNT);
#endif
}

/**
 *  Create data descriptor
 */
void HICMA_RUNTIME_desc_create(HICMA_desc_t *desc) {
    int64_t lmt = desc->lmt;
    int64_t lnt = desc->lnt;

    desc->occurences = 1;

    /*
     * Allocate starpu_handle_t array (handlers are initialized on the fly when
     * discovered by any algorithm to save space)
     */
    desc->schedopt = (void *) calloc(lnt * lmt, sizeof(starpu_data_handle_t));
    assert(desc->schedopt);

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
    /*
     * Register allocated memory as CUDA pinned memory
     */
    if ( (desc->use_mat == 1) && (desc->register_mat == 1) )
    {
        int64_t eltsze = HICMA_Element_Size(desc->dtyp);
        size_t size = (size_t)(desc->llm) * (size_t)(desc->lln) * eltsze;
        cudaError_t rc;

        /* Register the matrix as pinned memory */
        rc = cudaHostRegister( desc->mat, size, cudaHostRegisterPortable );
        if ( rc != cudaSuccess )
        {
            /* Disable the unregister as register failed */
            desc->register_mat = 0;
            hicma_warning("HiCMA_cham_RUNTIME_desc_create(StarPU): cudaHostRegister - ", cudaGetErrorString( rc ));
        }
    }
#endif

    if (desc->ooc) {
        int lastmm = desc->lm - (desc->lmt - 1) * desc->mb;
        int lastnn = desc->ln - (desc->lnt - 1) * desc->nb;
        int64_t eltsze = HICMA_Element_Size(desc->dtyp);
        int pagesize = getpagesize();

        if (((desc->mb * desc->nb * eltsze) % pagesize != 0) ||
            ((lastmm * desc->nb * eltsze) % pagesize != 0) ||
            ((desc->mb * lastnn * eltsze) % pagesize != 0) ||
            ((lastmm * lastnn * eltsze) % pagesize != 0)) {
            hicma_error("HiCMA_cham_RUNTIME_desc_create",
                        "Matrix and tile size not suitable for out-of-core: all tiles have to be multiples of 4096. Tip : choose 'n' and 'nb' as both multiples of 32.");
            return;
        }
    }

/*#if defined(HICMA_USE_MPI)
    {
        chameleon_starpu_tag_init( tag_width, tag_sep );

        if ( ((uintptr_t)(lnt*lmt)) > ((uintptr_t)(1UL<<tag_sep)) ) {
            hicma_fatal_error("HiCMA_cham_RUNTIME_desc_create", "Too many tiles in the descriptor for MPI tags");
            return;
        }
        assert(lmt*lmt<=(1<<tag_sep));

        if ( ((uintptr_t)desc->id) >= (uintptr_t)(1UL<<(tag_width-tag_sep)) ) {
            hicma_fatal_error("HiCMA_cham_RUNTIME_desc_create", "Number of descriptor available in MPI mode out of stock");
            return;
        }
        assert( ((uintptr_t)desc->id) < (uintptr_t)(1UL<<(tag_width-tag_sep)) );
    }
#endif*/
}

/**
 *  Destroy data descriptor
 */
void HICMA_RUNTIME_desc_destroy(HICMA_desc_t *desc) {
    desc->occurences--;

    /*
     * If this is the last descriptor using the matrix, we release the handle
     * and unregister the GPU data
     */
    if (desc->occurences == 0) {
        starpu_data_handle_t *handle = (starpu_data_handle_t *) (desc->schedopt);
        int lmt = desc->lmt;
        int lnt = desc->lnt;
        int m, n;

        for (n = 0; n < lnt; n++) {
            for (m = 0; m < lmt; m++) {
                if (*handle == NULL) {
                    handle++;
                    continue;
                }
                starpu_data_unregister(*handle);
                handle++;
            }
        }

#if defined(CHAMELEON_USE_CUDA) && !defined(CHAMELEON_SIMULATION)
        if ( (desc->use_mat == 1) && (desc->register_mat == 1) )
        {
            /* Unmap the pinned memory associated to the matrix */
            if (cudaHostUnregister(desc->mat) != cudaSuccess)
            {
                hicma_warning("HiCMA_cham_RUNTIME_desc_destroy(StarPU)",
                              "cudaHostUnregister failed to unregister the "
                              "pinned memory associated to the matrix");
            }
        }
#endif /* defined(CHAMELEON_USE_CUDA) */

        free(desc->schedopt);
    }
}

/**
 *  Acquire data
 */
int HICMA_RUNTIME_desc_acquire(const HICMA_desc_t *desc) {
    starpu_data_handle_t *handle = (starpu_data_handle_t *) (desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++) {
        for (m = 0; m < lmt; m++) {
            if ((*handle == NULL) ||
                !hicma_desc_islocal(desc, m, n)) {
                handle++;
                continue;
            }
            starpu_data_acquire(*handle, STARPU_R);
            handle++;
        }
    }
    return HICMA_SUCCESS;
}

/**
 *  Release data
 */
int HICMA_RUNTIME_desc_release(const HICMA_desc_t *desc) {
    starpu_data_handle_t *handle = (starpu_data_handle_t *) (desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++) {
        for (m = 0; m < lmt; m++) {
            if ((*handle == NULL) ||
                !hicma_desc_islocal(desc, m, n)) {
                handle++;
                continue;
            }
            starpu_data_release(*handle);
            handle++;
        }
    }
    return HICMA_SUCCESS;
}

/**
 *  Flush cached data
 */
void HICMA_RUNTIME_flush() {
#if defined(HICMA_USE_MPI)
    starpu_mpi_cache_flush_all_data(MPI_COMM_WORLD);
#endif
}

/**
 * Different implementations of the flush call based on StarPU version
 */
#ifdef HAVE_STARPU_DATA_WONT_USE

static inline void
hicma_chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_wont_use( handle );
}

#elif defined HAVE_STARPU_IDLE_PREFETCH

static inline void
chameleon_starpu_data_flush( starpu_data_handle_t handle)
{
    starpu_data_idle_prefetch_on_node(handle, STARPU_MAIN_RAM, 1);
    starpu_data_release_on_node(handle, -1);
}

static inline void
hicma_chameleon_starpu_data_wont_use( starpu_data_handle_t handle ) {
    starpu_data_acquire_on_node_cb( handle, -1, STARPU_R,
                                    chameleon_starpu_data_flush, handle );
}

#else

static inline void
hicma_chameleon_starpu_data_wont_use(starpu_data_handle_t handle) {
    starpu_data_acquire_cb(handle, STARPU_R,
                           (void (*)(void *)) &starpu_data_release, handle);
}

#endif

void HICMA_RUNTIME_desc_flush(const HICMA_desc_t *desc, const HICMA_sequence_t *sequence) {
    starpu_data_handle_t *handle = (starpu_data_handle_t *) (desc->schedopt);
    int lmt = desc->lmt;
    int lnt = desc->lnt;
    int m, n;

    for (n = 0; n < lnt; n++) {
        for (m = 0; m < lmt; m++, handle++) {
            if (*handle == NULL) {
                continue;
            }

#if defined(HICMA_USE_MPI)
            starpu_mpi_cache_flush( MPI_COMM_WORLD, *handle );
#endif
            if (hicma_desc_islocal(desc, m, n)) {
                hicma_chameleon_starpu_data_wont_use(*handle);
            }
        }
    }

    (void) sequence;
}

void HICMA_RUNTIME_data_flush(const HICMA_sequence_t *sequence,
                                   const HICMA_desc_t *A, int m, int n) {
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *handle = A->schedopt;
    handle += ((int64_t) A->lmt) * nn + mm;

    if (*handle == NULL) {
        return;
    }

#if defined(HICMA_USE_MPI)
    starpu_mpi_cache_flush( MPI_COMM_WORLD, *handle );
#endif

    if (hicma_desc_islocal(A, m, n)) {
        hicma_chameleon_starpu_data_wont_use(*handle);
    }

    (void) sequence;
}

#if defined(CHAMELEON_USE_MIGRATE)
void HICMA_RUNTIME_data_migrate( const HICMA_sequence_t *sequence,
                           const HICMA_desc_t *A, int Am, int An, int new_rank )
{
#if defined(HAVE_STARPU_MPI_DATA_MIGRATE)
    starpu_data_handle_t *handle = (starpu_data_handle_t*)(A->schedopt);
    starpu_data_handle_t lhandle;
    handle += ((int64_t)(A->lmt) * (int64_t)An + (int64_t)Am);

    lhandle = *handle;
    if ( lhandle == NULL ) {
        /* Register the data */
        lhandle = HICMA_RUNTIME_data_getaddr( A, Am, An );
    }

    starpu_mpi_data_migrate( MPI_COMM_WORLD, lhandle, new_rank );

    (void)sequence;
#else
    (void)sequence; (void)A; (void)Am; (void)An; (void)new_rank;
#endif
}
#endif

/**
 *  Get data addr
 */
/* For older revision of StarPU, STARPU_MAIN_RAM is not defined */
#ifndef STARPU_MAIN_RAM
#define STARPU_MAIN_RAM 0
#endif

void *HICMA_RUNTIME_data_getaddr(const HICMA_desc_t *A, int m, int n) {
    int64_t mm = m + (A->i / A->mb);
    int64_t nn = n + (A->j / A->nb);

    starpu_data_handle_t *ptrtile = A->schedopt;
    ptrtile += ((int64_t) A->lmt) * nn + mm;

    if (*ptrtile == NULL) {
        int home_node = -1;
        void *user_ptr = NULL;
        int myrank = A->myrank;
        int owner = A->get_rankof(A, m, n);
        int64_t eltsze = HICMA_Element_Size(A->dtyp);
        int tempmm = (mm == A->lmt - 1) ? (A->lm - mm * A->mb) : A->mb;
        int tempnn = (nn == A->lnt - 1) ? (A->ln - nn * A->nb) : A->nb;

        if (myrank == owner) {
            user_ptr = A->get_blkaddr(A, m, n);
            if (user_ptr != NULL) {
                home_node = STARPU_MAIN_RAM;
            }
        }

        starpu_matrix_data_register(ptrtile, home_node, (uintptr_t) user_ptr,
                                    BLKLDD(A, m),
                                    tempmm, tempnn, eltsze);

#ifdef HAVE_STARPU_DATA_SET_COORDINATES
        starpu_data_set_coordinates( *ptrtile, 2, m, n );
#endif

#if defined(HICMA_USE_MPI)
        {
            int64_t block_ind = A->lmt * nn + mm;
            starpu_mpi_data_register(*ptrtile, (A->id << tag_sep) | (block_ind), owner);
        }
#endif /* defined(HICMA_USE_MPI) */
    }

    return *ptrtile;
}
