/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 */

#include "hicma_struct.h"
static struct hicma_context hicma_context = {
    0, '\0', 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
};
const static struct hicma_context hicma_context_default = {
    0, '\0', 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
};
#include "flop_util_structs.h"
flop_counter counters[FLOP_NUMTHREADS];
