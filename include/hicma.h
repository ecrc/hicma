/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */
/**
 *
 * @file hicma.h
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon main header
 *
 * @version 1.0.0
 * @author Mathieu Faverge
 * @author Cedric Augonnet
 * @author Cedric Castagnede
 * @author Florent Pruvost
 * @date 2012-09-15
 *
 */


#ifndef __HICMA__
#define __HICMA__

#include <hicma_config.h>
#include <hicma_types.h>
#include <hicma_struct.h>
#include "hicma_constants.h"
#include "hicma_kernels.h"
#include "hicma_runtime.h"
#include "starsh.h"

int HICMA_get_print_mat();
int HICMA_set_print_mat();
int HICMA_get_print_index();
int HICMA_get_print_index_end();
int HICMA_set_print_index();
int HICMA_set_print_index_end();
int HICMA_unset_print_index_end();
int HICMA_set_use_fast_hcore_gemm();
int HICMA_get_use_fast_hcore_gemm();
int HICMA_get_always_fixed_rank();
int HICMA_get_fixed_rank();
int HICMA_set_fixed_rank(int rank);
int HICMA_set_starsh_format(STARSH_blrf *starsh_format);
STARSH_blrf* HICMA_get_starsh_format();

/** Morse header files */
/* Auxiliary */
int HICMA_Version           (int *ver_major, int *ver_minor, int *ver_micro);
int HICMA_My_Mpi_Rank       (void);
int HICMA_Init              (int nworkers, int ncudas);
int HICMA_InitPar           (int nworkers, int ncudas, int nthreads_per_worker);
int HICMA_Finalize          (void);
int HICMA_Pause             (void);
int HICMA_Resume            (void);
int HICMA_Distributed_start (void);
int HICMA_Distributed_stop  (void);
int HICMA_Comm_size         (void);
int HICMA_Comm_rank         (void);
int HICMA_Lapack_to_Tile    (void *Af77, int LDA, HICMA_desc_t *A);
int HICMA_Tile_to_Lapack    (HICMA_desc_t *A, void *Af77, int LDA);
int HICMA_Distributed_start (void);
int HICMA_Distributed_stop  (void);
int HICMA_Distributed_size  (int *size);
int HICMA_Distributed_rank  (int *rank);
int HICMA_GetThreadNbr      (void);

/* Descriptor */
int HICMA_Element_Size(int type);
int HICMA_Desc_Create  (HICMA_desc_t **desc, void *mat, HICMA_enum dtyp,
                        int mb, int nb, int bsiz, int lm, int ln,
                        int i, int j, int m, int n, int p, int q);
int HICMA_Desc_Create_User(HICMA_desc_t **desc, void *mat, HICMA_enum dtyp, int mb, int nb, int bsiz,
                           int lm, int ln, int i, int j, int m, int n, int p, int q,
                           void* (*get_blkaddr)( const HICMA_desc_t*, int, int ),
                           int (*get_blkldd)( const HICMA_desc_t*, int ),
                           int (*get_rankof)( const HICMA_desc_t*, int, int ));
int HICMA_Desc_Create_OOC(HICMA_desc_t **desc, HICMA_enum dtyp,
                          int mb, int nb, int bsiz, int lm, int ln,
                          int i, int j, int m, int n, int p, int q);
int HICMA_Desc_Create_OOC_User(HICMA_desc_t **desc, HICMA_enum dtyp,
                               int mb, int nb, int bsiz, int lm, int ln,
                               int i, int j, int m, int n, int p, int q,
                               int (*get_rankof)( const HICMA_desc_t*, int, int ));
int HICMA_Desc_Destroy (HICMA_desc_t **desc);
int HICMA_Desc_Acquire (HICMA_desc_t  *desc);
int HICMA_Desc_Release (HICMA_desc_t  *desc);
int HICMA_Desc_Flush   (HICMA_desc_t  *desc, HICMA_sequence_t *sequence);
void HICMA_user_tag_size(int, int) ;

/* Workspaces */
int HICMA_Dealloc_Workspace (HICMA_desc_t **desc);

/* Options */
int HICMA_Enable  (HICMA_enum option);
int HICMA_Disable (HICMA_enum option);
int HICMA_Set     (HICMA_enum param, int  value);
int HICMA_Get     (HICMA_enum param, int *value);
int HICMA_Set_HICMA_update_progress_callback(void (*p)(int, int)) ;

/* Sequences */
int HICMA_Sequence_Create  (HICMA_sequence_t **sequence);
int HICMA_Sequence_Destroy (HICMA_sequence_t *sequence);
int HICMA_Sequence_Wait    (HICMA_sequence_t *sequence);

#endif
