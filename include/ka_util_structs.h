#ifndef __KA_UTIL_STRUCTS__
#define __KA_UTIL_STRUCTS__
#include <stdint.h>
//printf("%s %s:%s:%d\t\t", ts[me],__FILE__, __FUNCTION__, __LINE__); \

typedef struct ka_time_t {
  unsigned long int __pad1;
  unsigned long int __pad2;
  int numthreads;
  int tid;
  double chol;
  double trsm;
  double syrk;
  double gemm;
  unsigned long int __pad3;
  unsigned long int __pad4;
} ka_time;
typedef struct ka_counter_t {
  unsigned long int __pad1;
  unsigned long int __pad2;
  unsigned long int copy;
  unsigned long int flop;
  int tid;
  unsigned long int potrf;
  unsigned long int trsm;
  unsigned long int syrk;
  unsigned long int update;
  unsigned long int geqrfu;
  unsigned long int geqrfv;
  unsigned long int geqrfv_gemm;
  unsigned long int gesvd;
  unsigned long int gesvd_gemm;
  unsigned long int gesvd_lasetu;
  unsigned long int gesvd_lasetv;
  unsigned long int gesvd_lacpyu;
  unsigned long int gesvd_lacpyv;
  unsigned long int gesvd_dscal;
  unsigned long int ormqru;
  unsigned long int ormqru_laset;
  unsigned long int ormqru_lacpy;
  unsigned long int ormqrv;
  unsigned long int ormqrv_laset;
  unsigned long int ormqrv_lacpy;
  unsigned long int omatcopy;

  int numthreads;
  unsigned long int __pad3;
  unsigned long int __pad4;
} ka_counter;
#define KA_NUMTHREADS 256

#ifdef KA_PRINT
#define ka_print(...) \
  do {\
  fflush(stdout);\
  fflush(stderr);\
  int me=omp_get_thread_num();\
  int nt=omp_get_num_threads();\
  int y;\
  ts[me][0]='\0';\
  for(y=0;y<me;y++)\
    strcat(ts[me], "\t");\
  sprintf(ts[me], "%s%d",ts[me],me);\
  for(y=0;y<nt-me;y++)\
    strcat(ts[me], "\t");\
  printf("%s\t", ts[me]); \
  printf(__VA_ARGS__); \
  fflush(stdout);\
  } while (0);
#else
#define ka_print(...) 
#endif
#if defined KA_DEBUG || defined KA_DEBUG_PRINTF
#define ka_printf(...) \
  do {\
  fflush(stdout);\
  fflush(stderr);\
  printf("%d %s:%s:%d\t\t", omp_get_thread_num(),__FILE__, __FUNCTION__, __LINE__); \
  printf(__VA_ARGS__); \
  fflush(stdout);\
  } while (0);
#else
#define ka_printf(...) \
  do{}while(0);  
#endif

#define ka_mark \
  fflush(stdout);\
  fflush(stderr);\
  printf("!!!%s:%d!!!",__FILE__, __LINE__);

#define ka_ex \
  fflush(stdout);\
  fflush(stderr);\
  printf("EARLY EXIT %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); \
  exit(0);

#ifdef KA_PAUSE
#define ka_pa \
  fflush(stdout);\
  fflush(stderr);\
  printf("PAUSED %s:%s:%d\n", __FILE__, __FUNCTION__, __LINE__); \
  getc(stdin);
#else
#define ka_pa \
 do{}while(0); 
#endif



#endif
