#ifndef __HICMA__
#define __HICMA__
#include "hicma_struct.h"
int HICMA_init();
int HICMA_get_print_mat();
int HICMA_set_print_mat();
int HICMA_get_print_index();
int HICMA_get_print_index_end();
int HICMA_set_print_index();
int HICMA_set_print_index_end();
int HICMA_unset_print_index_end();
int HICMA_set_use_fast_hcore_zgemm();
int HICMA_set_starsh_format(STARSH_blrf *starsh_format);
STARSH_blrf* HICMA_get_starsh_format();
int HICMA_get_use_fast_hcore_zgemm();
int HICMA_get_always_fixed_rank();
int HICMA_get_fixed_rank();
int HICMA_set_fixed_rank(int rank);
#endif
