/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file descriptor.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptors routines
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2010-11-15
 *
 ***
 *
 * @defgroup Descriptor
 * @brief Group descriptor routines exposed to users
 *
 */
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "hicma_descriptor.h"
#include <hicma_runtime.h>

static int nbdesc = 0;

/**
 ******************************************************************************
 *
 * @ingroup Descriptor
 *
 * hicma_desc_init_user - Internal function to create tiled matrix descriptor
 * with generic function for data distribution and storage format.
 *
 ******************************************************************************
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg HicmaRealFloat:     single precision real (S),
 *          @arg HicmaRealDouble:    double precision real (D),
 *          @arg HicmaComplexFloat:  single precision complex (C),
 *          @arg HicmaComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in bytes including padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 * @param[in] (*get_blkaddr)( const HICMA_desc_t *A, int m, int n)
 *          A function which return the address of the data corresponding to
 *          the tile A(m,n).
 *
 * @param[in] (*get_blkldd)( const HICMA_desc_t *A, int m )
 *          A function that return the leading dimension of the tile A(m,*).
 *
 * @param[in] (*get_rankof)( const HICMA_desc_t *A, int m, int n)
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @return  The descriptor with the matrix description parameters set.
 *
 */
HICMA_desc_t hicma_desc_init_user(HICMA_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j,
                                  int m,  int n,  int p, int q,
                                  void* (*get_blkaddr)( const HICMA_desc_t*, int, int ),
                                  int   (*get_blkldd) ( const HICMA_desc_t*, int      ),
                                  int   (*get_rankof) ( const HICMA_desc_t*, int, int ))
{
    HICMA_context_t *hicma;
    HICMA_desc_t desc;

    memset( &desc, 0, sizeof(HICMA_desc_t) );

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Desc_Create", "HiCMA not initialized");
        return desc;
    }

    // If one of the function get_* is NULL, we switch back to the default, like in hicma_desc_init()
    desc.get_blkaddr = get_blkaddr ? get_blkaddr : hicma_getaddr_ccrb;
    desc.get_blkldd  = get_blkldd  ? get_blkldd  : hicma_getblkldd_ccrb;
    desc.get_rankof  = get_rankof  ? get_rankof  : hicma_getrankof_2d;
    // Matrix properties
    desc.dtyp = dtyp;
    // Should be given as parameter to follow get_blkaddr (unused)
    desc.styp = HicmaCCRB;
    desc.mb   = mb;
    desc.nb   = nb;
    desc.bsiz = bsiz;
    // Large matrix parameters
    desc.lm = lm;
    desc.ln = ln;
    // Large matrix derived parameters
    desc.lmt = (lm%mb==0) ? (lm/mb) : (lm/mb+1);
    desc.lnt = (ln%nb==0) ? (ln/nb) : (ln/nb+1);
    // Submatrix parameters
    desc.i = i;
    desc.j = j;
    desc.m = m;
    desc.n = n;
    // Submatrix derived parameters
    desc.mt = (m == 0) ? 0 : (i+m-1)/mb - i/mb + 1;
    desc.nt = (n == 0) ? 0 : (j+n-1)/nb - j/nb + 1;

    desc.id = nbdesc; nbdesc++;
    desc.occurences = 0;
    desc.use_mat      = 1;
    desc.alloc_mat    = 1;
    desc.register_mat = (hicma->ncudas > 0) ? 1 : 0;
    desc.ooc          = 0;

    desc.myrank = HICMA_RUNTIME_comm_rank( hicma );

    // Grid size
    desc.p = p;
    desc.q = q;

    // Local dimensions in tiles
    if ( desc.myrank < (p*q) ) {
        desc.llmt = (desc.lmt + p - 1) / p;
        desc.llnt = (desc.lnt + q - 1) / q;

        // Local dimensions
        if ( ((desc.lmt-1) % p) == (desc.myrank / q) ) {
            desc.llm  = ( desc.llmt - 1 ) * mb + ((lm%mb==0) ? mb : (lm%mb));
        } else {
            desc.llm  =  desc.llmt * mb;
        }

        if ( ((desc.lnt-1) % q) == (desc.myrank % q) ) {
            desc.lln  = ( desc.llnt - 1 ) * nb + ((ln%nb==0) ? nb : (ln%nb));
        } else {
            desc.lln  =  desc.llnt * nb;
        }

        desc.llm1 = (desc.llm/mb);
        desc.lln1 = (desc.lln/nb);
    } else {
        desc.llmt = 0;
        desc.llnt = 0;
        desc.llm  = 0;
        desc.lln  = 0;
        desc.llm1 = 0;
        desc.lln1 = 0;
    }

    // Matrix address
    desc.mat = NULL;
    desc.A21 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb);
    desc.A12 = (size_t)(           desc.llm%mb)*(size_t)(desc.lln - desc.lln%nb) + desc.A21;
    desc.A22 = (size_t)(desc.llm - desc.llm%mb)*(size_t)(           desc.lln%nb) + desc.A12;

    return desc;
}

/**
 *  Internal static descriptor initializer
 */
HICMA_desc_t hicma_desc_init(HICMA_enum dtyp, int mb, int nb, int bsiz,
                             int lm, int ln, int i, int j,
                             int m,  int n,  int p, int q)
{
    return hicma_desc_init_user(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q,
                                hicma_getaddr_ccrb, hicma_getblkldd_ccrb, hicma_getrankof_2d);
}

/**
 *  Internal static descriptor initializer for a block diagonal matrix
 */
HICMA_desc_t hicma_desc_init_diag(HICMA_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j,
                                  int m,  int n,  int p, int q)
{
    return hicma_desc_init_user(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q,
                                hicma_getaddr_ccrb, hicma_getblkldd_ccrb, hicma_getrankof_2d_diag);
}

/**
 *  Internal static descriptor initializer for submatrices
 */
HICMA_desc_t* hicma_desc_submatrix(HICMA_desc_t *descA, int i, int j, int m, int n)
{
    HICMA_desc_t *descB = malloc(sizeof(HICMA_desc_t));
    int mb, nb;

    if ( (descA->i + i + m) > descA->lm ) {
        hicma_error("hicma_desc_submatrix", "The number of rows (i+m) of the submatrix doesn't fit in the parent matrix");
        assert((descA->i + i + m) > descA->lm);
    }
    if ( (descA->j + j + n) > descA->ln ) {
        hicma_error("hicma_desc_submatrix", "The number of rows (j+n) of the submatrix doesn't fit in the parent matrix");
        assert((descA->j + j + n) > descA->ln);
    }

    memcpy( descB, descA, sizeof(HICMA_desc_t) );
    mb = descA->mb;
    nb = descA->nb;
    // Submatrix parameters
    descB->i = descA->i + i;
    descB->j = descA->j + j;
    descB->m = m;
    descB->n = n;
    // Submatrix derived parameters
    descB->mt = (m == 0) ? 0 : (descB->i+m-1)/mb - descB->i/mb + 1;
    descB->nt = (n == 0) ? 0 : (descB->j+n-1)/nb - descB->j/nb + 1;

    // Increase the number of occurences to avoid multiple free of runtime specific data structures.
    descB->occurences++;

    return descB;
}

/**
 *  Check for descriptor correctness
 */
int hicma_desc_check(const HICMA_desc_t *desc)
{
    if (desc == NULL) {
        hicma_error("hicma_desc_check", "NULL descriptor");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    if (desc->mat == NULL && desc->use_mat == 1) {
        hicma_error("hicma_desc_check", "NULL matrix pointer");
        return HICMA_ERR_UNALLOCATED;
    }
    if (desc->dtyp != HicmaRealFloat &&
        desc->dtyp != HicmaRealDouble &&
        desc->dtyp != HicmaComplexFloat &&
        desc->dtyp != HicmaComplexDouble  ) {
        hicma_error("hicma_desc_check", "invalid matrix type");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->mb <= 0 || desc->nb <= 0) {
        hicma_error("hicma_desc_check", "negative tile dimension");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->bsiz < desc->mb*desc->nb) {
        hicma_error("hicma_desc_check", "tile memory size smaller than the product of dimensions");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->lm <= 0 || desc->ln <= 0) {
        hicma_error("hicma_desc_check", "negative matrix dimension");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->lm < desc->m) || (desc->ln < desc->n)) {
        hicma_error("hicma_desc_check", "matrix dimensions larger than leading dimensions");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if ((desc->i > 0 && desc->i >= desc->lm) || (desc->j > 0 && desc->j >= desc->ln)) {
        hicma_error("hicma_desc_check", "beginning of the matrix out of scope");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    if (desc->i+desc->m > desc->lm || desc->j+desc->n > desc->ln) {
        hicma_error("hicma_desc_check", "submatrix out of scope");
        return HICMA_ERR_ILLEGAL_VALUE;
    }
    return HICMA_SUCCESS;
}

/**
 *
 */
int hicma_desc_mat_alloc( HICMA_desc_t *desc )
{
    size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
        * (size_t)HICMA_Element_Size(desc->dtyp);
    if ((desc->mat = HICMA_RUNTIME_malloc(size)) == NULL) {
        hicma_error("hicma_desc_mat_alloc", "malloc() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }

    /* The matrix has already been registered by the Runtime alloc */
    desc->register_mat = 0;

    return HICMA_SUCCESS;
}

/**
 *
 */
int hicma_desc_mat_free( HICMA_desc_t *desc )
{
    if ( (desc->mat       != NULL) &&
         (desc->use_mat   == 1   ) &&
         (desc->alloc_mat == 1   ) )
    {
        size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
            * (size_t)HICMA_Element_Size(desc->dtyp);

        HICMA_RUNTIME_free(desc->mat, size);
        desc->mat = NULL;
    }
    return HICMA_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Create - Create tiled matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg HicmaRealFloat:     single precision real (S),
 *          @arg HicmaRealDouble:    double precision real (D),
 *          @arg HicmaComplexFloat:  single precision complex (C),
 *          @arg HicmaComplexDouble: double precision complex (Z).
 *
 * @param[in] mb
 *          Number of rows in a tile.
 *
 * @param[in] nb
 *          Number of columns in a tile.
 *
 * @param[in] bsiz
 *          Size in bytes including padding.
 *
 * @param[in] lm
 *          Number of rows of the entire matrix.
 *
 * @param[in] ln
 *          Number of columns of the entire matrix.
 *
 * @param[in] i
 *          Row index to the beginning of the submatrix.
 *
 * @param[in] j
 *          Column indes to the beginning of the submatrix.
 *
 * @param[in] m
 *          Number of rows of the submatrix.
 *
 * @param[in] n
 *          Number of columns of the submatrix.
 *
 * @param[in] p
 *          2D-block cyclic distribution in rows.
 *
 * @param[in] q
 *          2D-block cyclic distribution in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Create(HICMA_desc_t **descptr, void *mat, HICMA_enum dtyp, int mb, int nb, int bsiz,
                      int lm, int ln, int i, int j, int m, int n, int p, int q)
{
    HICMA_context_t *hicma;
    HICMA_desc_t *desc;
    int status;

    *descptr = NULL;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Desc_Create", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the descriptor */
    desc = (HICMA_desc_t*)malloc(sizeof(HICMA_desc_t));
    if (desc == NULL) {
        hicma_error("HiCMA_Desc_Create", "malloc() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }
    *desc = hicma_desc_init(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q);

    if (mat == NULL) {

        size_t size = (size_t)(desc->llm) * (size_t)(desc->lln)
            * (size_t)HICMA_Element_Size(desc->dtyp);

        if ((desc->mat = HICMA_RUNTIME_malloc(size)) == NULL) {
            hicma_error("HiCMA_Desc_Create", "malloc() failed");
            free(desc);
            return HICMA_ERR_OUT_OF_RESOURCES;
        }
        desc->use_mat      = 1;
        desc->alloc_mat    = 1;
        desc->register_mat = 0;

    } else {
        desc->mat = mat;
        /* memory of the matrix is handled by users */
        desc->alloc_mat    = 0;
        desc->use_mat      = 1;
        desc->register_mat = 0;
    }

    /* Create scheduler structure like registering data */
    HICMA_RUNTIME_desc_create( desc );

    status = hicma_desc_check( desc );
    if (status != HICMA_SUCCESS) {
        hicma_error("HiCMA_Desc_Create", "invalid descriptor");
        HICMA_Desc_Destroy( &desc );
        return status;
    }

    *descptr = desc;
    return HICMA_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Create_User - Create generic tiled matrix descriptor for general
 *  applications.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] mat
 *          Memory location of the matrix. If mat is NULL, the space to store
 *          the data is automatically allocated by the call to the function.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg HicmaRealFloat:     single precision real (S),
 *          @arg HicmaRealDouble:    double precision real (D),
 *          @arg HicmaComplexFloat:  single precision complex (C),
 *          @arg HicmaComplexDouble: double precision complex (Z).
 *
 * @param[in] nb
 *          Number of rows and columns in a tile.
 *
 * @param[in] m
 *          Number of rows of the entire matrix.
 *
 * @param[in] n
 *          Number of columns of the entire matrix.
 *
 * @param[in] p
 *          2d-block cyclic partitioning, number of tiles in rows.
 *
 * @param[in] q
 *          2d-block cyclic partitioning, number of tiles in columns.
 *
 * @param[in] (*get_blkaddr)( const HICMA_desc_t *A, int m, int n)
 *          A function which return the address of the data corresponding to
 *          the tile A(m,n).
 *
 * @param[in] (*get_blkldd)( const HICMA_desc_t *A, int m)
 *          A function that return the leading dimension of the tile A(m,*).
 *
 * @param[in] (*get_rankof)( const HICMA_desc_t *A, int m, int n)
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Create_User(HICMA_desc_t **descptr, void *mat, HICMA_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const HICMA_desc_t*, int, int ),
                           int   (*get_blkldd) ( const HICMA_desc_t*, int      ),
                           int   (*get_rankof) ( const HICMA_desc_t*, int, int ))
{
    HICMA_context_t *hicma;
    HICMA_desc_t *desc;
    int status;

    *descptr = NULL;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Desc_Create_User", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }

    /* Allocate memory and initialize the descriptor */
    desc = (HICMA_desc_t*)malloc(sizeof(HICMA_desc_t));
    if (desc == NULL) {
        hicma_error("HiCMA_Desc_Create_User", "malloc() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }

    *desc = hicma_desc_init_user(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q,
                                 get_blkaddr, get_blkldd, get_rankof);

    /* if the user gives a pointer to the overall data (tiles) we can use it */
    desc->use_mat = (mat == NULL) ? 0 : 1;

    /* memory of the matrix is handled by the user */
    desc->alloc_mat = 0;

    /* users data can have multiple forms: let him register tiles */
    desc->register_mat = 0;

    desc->mat = mat;

    /* Create runtime specific structure like registering data */
    HICMA_RUNTIME_desc_create( desc );

    status = hicma_desc_check( desc );
    if (status != HICMA_SUCCESS) {
        hicma_error("HiCMA_Desc_Create_User", "invalid descriptor");
        HICMA_Desc_Destroy( &desc );
        return status;
    }

    *descptr = desc;
    return HICMA_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Create_OOC_User - Create matrix descriptor for tiled matrix which
 *  may not fit memory.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg HicmaRealFloat:     single precision real (S),
 *          @arg HicmaRealDouble:    double precision real (D),
 *          @arg HicmaComplexFloat:  single precision complex (C),
 *          @arg HicmaComplexDouble: double precision complex (Z).
 *
 * @param[in] nb
 *          Number of rows and columns in a tile.
 *
 * @param[in] m
 *          Number of rows of the entire matrix.
 *
 * @param[in] n
 *          Number of columns of the entire matrix.
 *
 * @param[in] p
 *          2d-block cyclic partitioning, number of tiles in rows.
 *
 * @param[in] q
 *          2d-block cyclic partitioning, number of tiles in columns.
 *
 * @param[in] (*get_rankof)( const HICMA_desc_t *A, int m, int n)
 *          A function that return the MPI rank of the tile A(m,n).
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Create_OOC_User(HICMA_desc_t **descptr, HICMA_enum dtyp, int mb, int nb, int bsiz,
                               int lm, int ln, int i, int j, int m, int n, int p, int q,
                               int (*get_rankof)( const HICMA_desc_t*, int, int ))
{
#if !defined (CHAMELEON_SCHED_STARPU)
    (void)descptr; (void)dtyp; (void)mb; (void)nb; (void)bsiz;
    (void)lm; (void)ln; (void)i; (void)j; (void)m; (void)n; (void)p; (void)q;
    (void)get_rankof;

    hicma_error("HiCMA_Desc_Create_OOC_User", "Only StarPU supports on-demand tile allocation");
    return HICMA_ERR_NOT_INITIALIZED;
#else
    HICMA_context_t *hicma;
    HICMA_desc_t *desc;
    int status;

    *descptr = NULL;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Desc_Create_OOC_User", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }
    /* Allocate memory and initialize the descriptor */
    desc = (HICMA_desc_t*)malloc(sizeof(HICMA_desc_t));
    if (desc == NULL) {
        hicma_error("HiCMA_Desc_Create_OOC_User", "malloc() failed");
        return HICMA_ERR_OUT_OF_RESOURCES;
    }
    *desc = hicma_desc_init_user(dtyp, mb, nb, bsiz, lm, ln, i, j, m, n, p, q,
                                 hicma_getaddr_null, NULL, get_rankof);

    /* memory of the matrix is completely handled by runtime */
    desc->use_mat      = 0;
    desc->alloc_mat    = 0;
    desc->register_mat = 0;

    desc->mat = NULL;
    desc->ooc = 1;

    /* Create scheduler structure like registering data */
    HICMA_RUNTIME_desc_create( desc );

    status = hicma_desc_check( desc );
    if (status != HICMA_SUCCESS) {
        hicma_error("HiCMA_Desc_Create_OOC_User", "invalid descriptor");
        HICMA_Desc_Destroy( &desc );
        return status;
    }

    *descptr = desc;
    return HICMA_SUCCESS;
#endif
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Create_OOC - Create matrix descriptor for tiled matrix which may
 *  not fit memory.
 *
 ******************************************************************************
 *
 * @param[out] desc
 *          On exit, descriptor of the matrix.
 *
 * @param[in] dtyp
 *          Data type of the matrix:
 *          @arg HicmaRealFloat:     single precision real (S),
 *          @arg HicmaRealDouble:    double precision real (D),
 *          @arg HicmaComplexFloat:  single precision complex (C),
 *          @arg HicmaComplexDouble: double precision complex (Z).
 *
 * @param[in] nb
 *          Number of rows and columns in a tile.
 *
 * @param[in] m
 *          Number of rows of the entire matrix.
 *
 * @param[in] n
 *          Number of columns of the entire matrix.
 *
 * @param[in] p
 *          2d-block cyclic partitioning, number of tiles in rows.
 *
 * @param[in] q
 *          2d-block cyclic partitioning, number of tiles in columns.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Create_OOC(HICMA_desc_t **descptr, HICMA_enum dtyp, int mb, int nb, int bsiz,
                          int lm, int ln, int i, int j, int m, int n, int p, int q)
{
    return HICMA_Desc_Create_OOC_User( descptr, dtyp, mb, nb, bsiz,
                                       lm, ln, i, j, m, n, p, q,
                                       hicma_getrankof_2d );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Destroy - Destroys matrix descriptor.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Destroy(HICMA_desc_t **desc)
{
    HICMA_context_t *hicma;

    hicma = hicma_context_self();
    if (hicma == NULL) {
        hicma_error("HiCMA_Desc_Destroy", "HiCMA not initialized");
        return HICMA_ERR_NOT_INITIALIZED;
    }

    if (*desc == NULL) {
        hicma_error("HiCMA_Desc_Destroy", "attempting to destroy a NULL descriptor");
        return HICMA_ERR_UNALLOCATED;
    }

    HICMA_RUNTIME_desc_destroy( *desc );
    hicma_desc_mat_free( *desc );
    free(*desc);
    *desc = NULL;
    return HICMA_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Acquire - Ensures that all data of the descriptor are
 *  up-to-date.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Acquire (HICMA_desc_t  *desc) {
    return HICMA_RUNTIME_desc_acquire( desc );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Release - Release the data of the descriptor acquired by the
 *  application. Should be called if HICMA_Desc_Acquire has been called on the
 *  descriptor and if you do not need to access to its data anymore.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Release (HICMA_desc_t  *desc) {
    return HICMA_RUNTIME_desc_release( desc );
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_Desc_Flush - Flushes the data in the sequence when they won't be reused. This calls cleans up the distributed communication caches, and transfer the data back to the CPU.
 *
 ******************************************************************************
 *
 * @param[in] desc
 *          Matrix descriptor.
 *
 ******************************************************************************
 *
 * @return
 *          \retval HICMA_SUCCESS successful exit
 *
 */
int HICMA_Desc_Flush( HICMA_desc_t     *desc,
                      HICMA_sequence_t *sequence )
{
    HICMA_RUNTIME_desc_flush( desc, sequence );
    return HICMA_SUCCESS;
}

/**
 *****************************************************************************
 *
 * @ingroup Descriptor
 *
 *  HICMA_user_tag_size - Set the sizes for the MPI tags
 *  Default value: tag_width=31, tag_sep=24, meaning that the MPI tag is stored
 *  in 31 bits, with 24 bits for the tile tag and 7 for the descriptor.  This
 *  function must be called before any descriptor creation.
 *
 ******************************************************************************
 *
 * @param[in] user_tag_width
 *          The new value for tag_width.
 *
 * @param[in] user_tag_sep
 *          The new value for tag_sep.
 *
 ******************************************************************************
 *
 * @return
 *          \retval none
 *
 */
void HICMA_user_tag_size(int user_tag_width, int user_tag_sep) {
    HICMA_RUNTIME_comm_set_tag_sizes( user_tag_width, user_tag_sep );
    return;
}
