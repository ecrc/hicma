#ifndef __KA_DEFINES__
#define __KA_DEFINES__
#define KA_PPROFILE 
//#define MAX_DIM 8 
//#define MAX_DIM 16 
//#define MAX_DIM 32 
#define MAX_DIM 8192 
//#define MAX_DIM 16384 
//#define MAX_DIM 32768 
//#define KA_DEBUG 
//#define KA_DEBUG_PRINTF
#define KA_PAUSE 
///Low rank block allocation strategy
#define KA_LW_EXACT 1
#define KA_LW_MORE  2
//#define KA_LW KA_LW_EXACT
#define KA_LW KA_LW_MORE
int KA_LW_SIZE_MORE_MEM = 0;
int print_ops = 0;
//#define KA_COUNT
#endif
