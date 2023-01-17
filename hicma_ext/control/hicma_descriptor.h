/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file descriptor.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon descriptor header
 *
 * @version 1.0.0
 * @author Jakub Kurzak
 * @author Mathieu Faverge
 * @author Cedric Castagnede
 * @date 2012-09-15
 *
 */
#ifndef _HICMA_CHAM_DESCRIPTOR_H_
#define _HICMA_CHAM_DESCRIPTOR_H_

#include <assert.h>
#include <hicma_config.h>
#include <hicma_struct.h>
#include <control/hicma_auxiliary.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 *  Internal routines
 */
inline static void* hicma_geteltaddr(const HICMA_desc_t *A, int m, int n, int eltsize);
inline static void* hicma_getaddr_cm    (const HICMA_desc_t *A, int m, int n);
inline static void* hicma_getaddr_ccrb  (const HICMA_desc_t *A, int m, int n);
inline static void* hicma_getaddr_null  (const HICMA_desc_t *A, int m, int n);
inline static int   hicma_getblkldd_cm  (const HICMA_desc_t *A, int m);
inline static int   hicma_getblkldd_ccrb(const HICMA_desc_t *A, int m);

/**
 *  Data distributions
 */
inline static int   hicma_getrankof_2d(const HICMA_desc_t *desc, int m, int n);
inline static int   hicma_getrankof_2d_diag(const HICMA_desc_t *desc, int m, int n);

HICMA_desc_t hicma_desc_init(HICMA_enum dtyp, int mb, int nb, int bsiz,
                             int lm, int ln, int i, int j, int m, int n, int p, int q);
HICMA_desc_t hicma_desc_init_diag(HICMA_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j, int m, int n, int p, int q);
HICMA_desc_t hicma_desc_init_user(HICMA_enum dtyp, int mb, int nb, int bsiz,
                                  int lm, int ln, int i, int j,
                                  int m,  int n,  int p, int q,
                                  void* (*get_blkaddr)( const HICMA_desc_t*, int, int ),
                                  int (*get_blkldd)( const HICMA_desc_t*, int ),
                                  int (*get_rankof)( const HICMA_desc_t*, int, int ));
HICMA_desc_t* hicma_desc_submatrix(HICMA_desc_t *descA, int i, int j, int m, int n);

int hicma_desc_check    (const HICMA_desc_t *desc);
int hicma_desc_mat_alloc(HICMA_desc_t *desc);
int hicma_desc_mat_free (HICMA_desc_t *desc);

#define BLKLDD(A, k) A->get_blkldd( A, k )

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static void* hicma_getaddr_ccrb(const HICMA_desc_t *A, int m, int n)
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = HICMA_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(HICMA_USE_MPI)
    //assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = (size_t)(A->bsiz) * (mm + (size_t)(A->llm1) * nn);
        else
            offset = A->A12 + ((size_t)(A->mb * (A->lln%A->nb)) * mm);
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((size_t)((A->llm%A->mb) * A->nb) * nn);
        else
            offset = A->A22;
    }

    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 */
inline static void *hicma_getaddr_cm(const HICMA_desc_t *A, int m, int n)
{
    size_t mm = m + A->i / A->mb;
    size_t nn = n + A->j / A->nb;
    size_t eltsize = HICMA_Element_Size(A->dtyp);
    size_t offset = 0;

#if defined(HICMA_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    offset = (size_t)(A->llm * A->nb) * nn + (size_t)(A->mb) * mm;
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 *  Internal function to return address of block (m,n) with m,n = block indices
 *  This version lets the runtime allocate on-demand.
 */
inline static void *hicma_getaddr_null(const HICMA_desc_t *A, int m, int n)
{
    (void)A; (void)m; (void)n;
    return NULL;
}

/**
 *  Internal function to return address of element A(m,n) with m,n = matrix indices
 */
inline static void* hicma_geteltaddr(const HICMA_desc_t *A, int m, int n, int eltsize) // Not used anywhere ?!
{
    size_t mm = (m + A->i)/A->mb;
    size_t nn = (n + A->j)/A->nb;
    size_t offset = 0;

#if defined(HICMA_USE_MPI)
    assert( A->myrank == A->get_rankof( A, mm, nn) );
    mm = mm / A->p;
    nn = nn / A->q;
#endif

    if (mm < (size_t)(A->llm1)) {
        if (nn < (size_t)(A->lln1))
            offset = A->bsiz*(mm+A->llm1*nn) + m%A->mb + A->mb*(n%A->nb);
        else
            offset = A->A12 + (A->mb*(A->lln%A->nb)*mm) + m%A->mb + A->mb*(n%A->nb);
    }
    else {
        if (nn < (size_t)(A->lln1))
            offset = A->A21 + ((A->llm%A->mb)*A->nb*nn) + m%A->mb + (A->llm%A->mb)*(n%A->nb);
        else
            offset = A->A22 + m%A->mb  + (A->llm%A->mb)*(n%A->nb);
    }
    return (void*)((intptr_t)A->mat + (offset*eltsize) );
}

/**
 *  Internal function to return the leading dimension of element A(m,*) with m,n = block indices
 */
inline static int hicma_getblkldd_ccrb(const HICMA_desc_t *A, int m)
{
    int mm = m + A->i / A->mb;
    return ( ((mm+1) == A->lmt) && ((A->lm % A->mb) != 0)) ? A->lm % A->mb : A->mb;
}

inline static int hicma_getblkldd_cm(const HICMA_desc_t *A, int m) {
    (void)m;
    return A->llm;
}


/**
 *  Internal function to return MPI rank of element A(m,n) with m,n = block indices
 */
inline static int hicma_getrankof_2d(const HICMA_desc_t *A, int m, int n)
{
    int mm = m + A->i / A->mb;
    int nn = n + A->j / A->nb;
    return (mm % A->p) * A->q + (nn % A->q);
}

/**
 *  Internal function to return MPI rank of element DIAG(m,0) with m,n = block indices
 */
inline static int hicma_getrankof_2d_diag(const HICMA_desc_t *A, int m, int n)
{
    int mm = m + A->i / A->mb;
    assert( n == 0 );
    return (mm % A->p) * A->q + (mm % A->q);
}


/**
 * Detect if the tile is local or not
 */
inline static int hicma_desc_islocal( const HICMA_desc_t *A, int m, int n )
{
#if defined(HICMA_USE_MPI)
    return (A->myrank == A->get_rankof(A, m, n));
#else
    (void)A; (void)m; (void)n;
    return 1;
#endif /* defined(HICMA_USE_MPI) */
}

/**
 * Declare data accesses of codelets using these macros, for instance:
 * HICMA_BEGIN_ACCESS_DECLARATION
 * HICMA_ACCESS_R(A, Am, An)
 * HICMA_ACCESS_R(B, Bm, Bn)
 * HICMA_ACCESS_RW(C, Cm, Cn)
 * HICMA_END_ACCESS_DECLARATION
 */
#define HICMA_BEGIN_ACCESS_DECLARATION { \
    unsigned __hicma_need_submit = 0; \
    RUNTIME_BEGIN_ACCESS_DECLARATION

#define HICMA_ACCESS_R(A, Am, An) do { \
    if (hicma_desc_islocal(A, Am, An)) __hicma_need_submit = 1; \
    RUNTIME_ACCESS_R(A, Am, An); \
} while(0)

#define HICMA_ACCESS_W(A, Am, An) do { \
    if (hicma_desc_islocal(A, Am, An)) __hicma_need_submit = 1; \
    RUNTIME_ACCESS_W(A, Am, An); \
} while(0)

#define HICMA_ACCESS_RW(A, Am, An) do { \
    if (hicma_desc_islocal(A, Am, An)) __hicma_need_submit = 1; \
    RUNTIME_ACCESS_RW(A, Am, An); \
} while(0)

#define HICMA_RANK_CHANGED(rank) do {\
    __hicma_need_submit = 1; \
    RUNTIME_RANK_CHANGED(rank); \
} while (0)

#define HICMA_END_ACCESS_DECLARATION \
    RUNTIME_END_ACCESS_DECLARATION; \
    if (!__hicma_need_submit) return; \
}

#ifdef __cplusplus
}
#endif

#endif
