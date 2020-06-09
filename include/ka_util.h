#ifndef __KA_UTIL__
#define __KA_UTIL__
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
void ka_util_init_counters(int numthreads){
  int i;
  for(i=0; i<KA_NUMTHREADS; i++){
    counters[i].copy = 0; 
    counters[i].flop = 0;
    counters[i].potrf = 0;
    counters[i].trsm = 0;
    counters[i].syrk = 0;
    counters[i].update = 0;
    counters[i].geqrfu = 0;
    counters[i].geqrfv = 0;
    counters[i].geqrfv_gemm = 0;
    counters[i].gesvd = 0;
    counters[i].gesvd_gemm = 0;
    counters[i].gesvd_lasetu = 0;
    counters[i].gesvd_lasetv = 0;
    counters[i].gesvd_lacpyu = 0;
    counters[i].gesvd_lacpyv = 0;
    counters[i].gesvd_dscal = 0;
    counters[i].ormqru = 0;
    counters[i].ormqru_laset = 0;
    counters[i].ormqru_lacpy = 0;
    counters[i].ormqrv = 0;
    counters[i].ormqrv_laset = 0;
    counters[i].ormqrv_lacpy = 0;
    counters[i].omatcopy = 0;
    counters[i].tid = i;
    counters[i].numthreads = numthreads;
  }
} 
ka_counter ka_util_sum_counters(int numthreads){
  ka_counter res;
  res.copy = 0;
  res.flop = 0;
  res.tid = 0;
  res.potrf         = 0;
  res.trsm          = 0;
  res.syrk          = 0;
  res.update        = 0;
  res.geqrfu        = 0;
  res.geqrfv        = 0;
  res.geqrfv_gemm   = 0;
  res.gesvd         = 0;
  res.gesvd_gemm    = 0;
  res.gesvd_lasetu  = 0;
  res.gesvd_lasetv  = 0;
  res.gesvd_lacpyu  = 0;
  res.gesvd_lacpyv  = 0;
  res.gesvd_dscal   = 0;
  res.ormqru        = 0;
  res.ormqru_laset  = 0;
  res.ormqru_lacpy  = 0;
  res.ormqrv        = 0;
  res.ormqrv_laset  = 0;
  res.ormqrv_lacpy  = 0;
  res.omatcopy      = 0;

  int i;
  for(i=0; i<numthreads; i++){
    res.copy         += counters[i].copy; 
    res.flop         += counters[i].flop;
    res.tid          += counters[i].tid;
    res.potrf        += counters[i].potrf;
    res.trsm         += counters[i].trsm;
    res.syrk         += counters[i].syrk;
    res.update       += counters[i].update;
    res.geqrfu       += counters[i].geqrfu;
    res.geqrfv       += counters[i].geqrfv;
    res.geqrfv_gemm  += counters[i].geqrfv_gemm;
    res.gesvd        += counters[i].gesvd;
    res.gesvd_gemm   += counters[i].gesvd_gemm;
    res.gesvd_lasetu += counters[i].gesvd_lasetu;
    res.gesvd_lasetv += counters[i].gesvd_lasetv;
    res.gesvd_lacpyu += counters[i].gesvd_lacpyu;
    res.gesvd_lacpyv += counters[i].gesvd_lacpyv;
    res.gesvd_dscal  += counters[i].gesvd_dscal;
    res.ormqru       += counters[i].ormqru;
    res.ormqru_laset += counters[i].ormqru_laset;
    res.ormqru_lacpy += counters[i].ormqru_lacpy;
    res.ormqrv       += counters[i].ormqrv;
    res.ormqrv_laset += counters[i].ormqrv_laset;
    res.ormqrv_lacpy += counters[i].ormqrv_lacpy;
    res.omatcopy     += counters[i].omatcopy;
  }
  return res;
} 
void ka_util_sprint_counter(ka_counter c, char* str){
  if(str == NULL){
    return; 
  }
  //sprintf(str, "%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu\t%lu"
  sprintf(str, "%lu\t%lu\t%lu\t%lu"
                                           ,c.potrf        
                                           ,c.trsm         
                                           ,c.syrk         
                                           ,c.update    /*    
                                           ,c.flop
                                           ,c.copy
                                           ,c.geqrfu       
                                           ,c.geqrfv       
                                           ,c.geqrfv_gemm  
                                           ,c.gesvd        
                                           ,c.gesvd_gemm   
                                           ,c.gesvd_lasetu 
                                           ,c.gesvd_lasetv 
                                           ,c.gesvd_lacpyu 
                                           ,c.gesvd_lacpyv 
                                           ,c.gesvd_dscal  
                                           ,c.ormqru       
                                           ,c.ormqru_lacpy 
                                           ,c.ormqru_lacpy 
                                           ,c.ormqrv       
                                           ,c.ormqrv_laset 
                                           ,c.ormqrv_lacpy 
                                           ,c.omatcopy*/    
                                           );
} 
void ka_util_print_counter(ka_counter c){
  printf( "copy:%lu\tflop:%lu\t%g gflops\tpotrf:%lu\ttrsm:%lu\tsyrk:%lu\tupdate:%lu\nqra:%lu\tomatcpy:%lu\tqrb:%lu\tqrbgemm:%lu\nsvd:%lu\tgemm:%lu\tsetu:%lu\tsetv:%lu\tcpu:%lu\tcpv:%lu\tdscal:%lu\normu:%lu\tset:%lu\tcp:%lu\normv:%lu\tset:%lu\tcp:%lu\n", c.copy
                                           ,c.flop  
                                           ,(c.flop/(1024.0*1024*1024)) //gflop  
                                           ,c.potrf        
                                           ,c.trsm         
                                           ,c.syrk         
                                           ,c.update         
                                           ,c.geqrfu       
                                           ,c.omatcopy    
                                           ,c.geqrfv       
                                           ,c.geqrfv_gemm  
                                           ,c.gesvd        
                                           ,c.gesvd_gemm   
                                           ,c.gesvd_lasetu 
                                           ,c.gesvd_lasetv 
                                           ,c.gesvd_lacpyu 
                                           ,c.gesvd_lacpyv 
                                           ,c.gesvd_dscal  
                                           ,c.ormqru       
                                           ,c.ormqru_laset 
                                           ,c.ormqru_lacpy 
                                           ,c.ormqrv       
                                           ,c.ormqrv_laset 
                                           ,c.ormqrv_lacpy 
                                           );
} 

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


void ka_util_print_array_rowmajor_double(int nelm, double* arr) {
  int i;
  for(i=0; i<nelm; i++){
    printf("%+.4f\t",arr[i]);
  }
  printf("\n");
}

uint32_t xor128(void) {
  static uint32_t x = 123456789;
  static uint32_t y = 362436069;
  static uint32_t z = 521288629;
  static uint32_t w = 88675123;
  uint32_t t;
  t = x ^ (x << 11);
  x = y; y = z; z = w;
  return w = w ^ (w >> 19) ^ (t ^ (t >> 8));
}
double ka_util_random_double(){
  //double X=((double)rand()/(double)RAND_MAX);
  double X=((double)xor128()/(double)RAND_MAX);
  return X;
}
void ka_util_print_2DarrRow(int m, int n, double* M){
  int i, j;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      printf("%.3g\t",M[i*m+j]);
    }
    printf("\n");
  }
}
void  ka_util_print_2DarrCol(int m, int n, double* M){
  int i, j;
  for(i = 0; i < m; i++){
    for(j = 0; j < n; j++){
      printf("%.3g\t",M[j*n+i]);
    }
    printf("\n");
  }
}
int ka_util_min2int(int a, int b){
  if(a<b){
    return a;
  } else {
    return b;
  }
}
int ka_util_max2int(int a, int b){
  if(a>b){
    return a;
  } else {
    return b;
  }
}

#endif
