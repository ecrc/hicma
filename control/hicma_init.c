#include "hicma_init.h"
int HICMA_init(){
    hicma_context = hicma_context_default;
    return 0;
}
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
int HICMA_set_use_fast_hcore_zgemm(){
    hicma_context.use_fast_hcore_zgemm = 1;
    return 0;
}
int HICMA_get_use_fast_hcore_zgemm(){
    return hicma_context.use_fast_hcore_zgemm;
}
int HICMA_unset_use_fast_hcore_zgemm(){
    hicma_context.use_fast_hcore_zgemm = 0;
    return 0;
}
int HICMA_set_starsh_format(STARSH_blrf *starsh_format){
    hicma_context.starsh_format = starsh_format; 
    return 0;
}
STARSH_blrf * HICMA_get_starsh_format(){
    return hicma_context.starsh_format; 
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
