#include "hicma_struct.h"
static struct hicma_context hicma_context = {
    0, '\0', 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
};
const static struct hicma_context hicma_context_default = {
    0, '\0', 0, 0 , 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 
};
//#define KA_COUNT
#ifdef KA_COUNT
#include "ka_defines.h"
#include "ka_util_structs.h"
ka_counter counters[KA_NUMTHREADS];
#endif
