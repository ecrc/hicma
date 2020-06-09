#ifndef __KA_COUNTS__
#define __KA_COUNTS__
unsigned long int ka_counts(char op, unsigned long int a, unsigned long int b, unsigned long int c, unsigned long int d){
  unsigned long int res = 0;
  if(op == 'q') {//geqrf  if m >= n
    unsigned long int m = a;
    unsigned long int n = b;
    res = 2*m*n*n - (unsigned long int)(2*n*n*n/3.0f) + 2*m*n + (unsigned long int)(17*n/3.0f);
    //printf("%lu %lu %lu\n", m, n, res);
  } 
  else if(op == 'c') {//potrf  
    unsigned long int n = a;
    res = n*n*n/3 - n*n/2.0 + n/6 ;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 't') {//trsm  
    unsigned long int m = a;
    unsigned long int n = b;
    int side = c; //1:left 2:right
    if(side == 1)
        res = n*m*m;
    else if(side == 2)
        res = m*n*n;
    else
        fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 'm') {//potrf  
    unsigned long int m = a;
    unsigned long int n = b;
    unsigned long int k = c;
    res = m*n*k*2;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 's') {//svd  
    unsigned long int n = a;
    res = 22*n*n*n;
    //printf("%lu %lu\n",  n, res);
  }
  else if(op == 'o') {//ormqr  
    unsigned long int m = a;
    unsigned long int n = b;
    unsigned long int k = c;
    int side = d; //1:left 2:right
    if(side == 1)
        res = 4*n*m*k-2*n*k*k+3*n*k;
    else if(side == 2)
        res = 4*n*m*k-2*m*k*k+2*m*k+n*k-k*k/2+k/2;
    else
        fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
    //printf("%lu %lu\n",  n, res);
  } else if (op == 'r') {//trmm
      /**
       * B is an m-by-n matrix,
       * A is a unit, or non-unit, upper or lower triangular matrix
       * left   B := alpha*op(A)*B
       * right  B := alpha*B*op(A)
       */
      unsigned long int m = a;
      unsigned long int n = b;
      int side = c;
      if(side == 1) //1:left 2:right
          res = n * m * m;
      else if (side == 2)
          res = m * n * n;
      else
          fprintf(stderr, "%s %d: invalid side:%d\n", __FILE__, __LINE__, side);
  }
  return res;
}

#endif
