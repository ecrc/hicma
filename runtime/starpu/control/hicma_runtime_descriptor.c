


#include <stdlib.h>
#include <unistd.h>
static int tag_sep   = 24;
void *HICMA_data_getaddr( const MORSE_desc_t *A,  int m, int n, int ncols)
{
	int64_t mm = m + (A->i / A->mb);
	int64_t nn = n + (A->j / A->nb);

	starpu_data_handle_t *ptrtile = A->schedopt;
	ptrtile += ((int64_t)A->lmt) * nn + mm;

	if (*ptrtile == NULL || phasenumber==2) {
		int home_node = -1;
		void *user_ptr = NULL;
		int myrank = A->myrank;
		int owner  = A->get_rankof( A, m, n );
		int64_t eltsze = MORSE_Element_Size(A->dtyp);
		int tempmm = (mm == A->lmt-1) ? (A->lm - mm * A->mb) : A->mb;
		//int tempnn = (nn == A->lnt-1) ? (A->ln - nn * A->nb) : A->nb;
                int tempnn = ncols;

//		printf(" %d -%d -%d\n", m, n, eltsze);

		if ( myrank == owner ) {
			user_ptr = A->get_blkaddr(A, m, n);
			if ( user_ptr != NULL ) {
				home_node = STARPU_MAIN_RAM;
			}
		}

		starpu_matrix_data_register( ptrtile, home_node, (uintptr_t) user_ptr,
				BLKLDD(A, m),
				tempmm, tempnn, eltsze );

#ifdef HAVE_STARPU_DATA_SET_COORDINATES
		starpu_data_set_coordinates( *ptrtile, 2, m, n );
#endif

#if defined(CHAMELEON_USE_MPI)
		{
			int64_t block_ind = A->lmt * nn + mm;
			starpu_mpi_data_register(*ptrtile, (A->id << tag_sep) | (block_ind), owner);
		}
#endif /* defined(MORSE_USE_MPI) */
	}

	return *ptrtile;
}
