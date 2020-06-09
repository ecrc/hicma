/**
 * @copyright (c) 2017 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_zgytlr.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.0
 * @author Kadir Akbudak
 * @date 2017-11-16
 * @precisions normal z -> c d s
 **/
#include "morse.h"
#include "runtime/starpu/chameleon_starpu.h"
//#include "runtime/starpu/include/runtime_codelet_z.h"
#include "hcore_z.h"
#include "hicma_common.h"
#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(gytlrb)
extern int isPrecentageused;
extern float precent1;
//CHAMELEON_CL_CB(zgytlr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N)
/*   MORSE_TASK_zgytlr - Generate a tile for random matrix. */

//#if defined(HICMA_COMPLEX)
void AL4SAN_Matrix_Generation(int *nip, int *ntrian, MORSE_Complex64_t *zz, int *q, int *p, int *local_nt, int *nb);

void HICMA_TASK_zgytlr_buffer( const MORSE_option_t *options,
                       int m, int n,
                       const MORSE_desc_t *AUV,
                       const MORSE_desc_t *Ark,
                       int Am, int An,
                       int lda,
                       int ldu,
                       int ldv,
                       int bigM, int m0, int n0, unsigned long long int seed,
                       int maxrank, double tol,
                       int compress_diag, int ntrian, int nip, int local_nt, int NB
//                       const MORSE_desc_t *Dense
                       )
{
    struct starpu_codelet *codelet = &cl_zgytlrb;
    //void (*callback)(void*) = options->profiling ? cl_zgytlr_callback : NULL;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;
    
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(AUV, Am, An);
//    MORSE_ACCESS_W(Dense, Am, An);
    MORSE_ACCESS_RW(Ark, Am, An);
    MORSE_END_ACCESS_DECLARATION;
    
    //printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);
    //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    //printf("\n nAUV:%d m:%d, n:%d\n", nAUV, m, n);
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,    &m,                      sizeof(int),
                       STARPU_VALUE,    &n,                      sizeof(int),
                       STARPU_VALUE,    &nAUV,                      sizeof(int),
                       STARPU_W,         RTBLKADDR(AUV, MORSE_Complex64_t, Am, An),
                       STARPU_RW,         RTBLKADDR(Ark, double, Am, An),
//                       STARPU_R,         RTBLKADDR(Dense, MORSE_Complex64_t, Am, An), // _R must be _W SERIOUSLY. BUT _W STALLS ON SHAHEEN. FIXME
                       STARPU_VALUE,  &lda,                      sizeof(int),
                       STARPU_VALUE,  &ldu,                      sizeof(int),
                       STARPU_VALUE,  &ldv,                      sizeof(int),
                       STARPU_VALUE, &bigM,                      sizeof(int),
                       STARPU_VALUE,   &Am,                      sizeof(int),
                       STARPU_VALUE,   &An,                      sizeof(int),
                       STARPU_VALUE, &seed,   sizeof(unsigned long long int),
                       STARPU_VALUE,   &maxrank,                      sizeof(int),
                       STARPU_VALUE,   &tol,                      sizeof(double),
                       STARPU_VALUE,   &compress_diag,            sizeof(int),
                       STARPU_VALUE,  &ntrian,                      sizeof(int),
                       STARPU_VALUE, &nip,                      sizeof(int),
                       STARPU_VALUE,   &local_nt,                      sizeof(int),
                       STARPU_VALUE,   &NB,                      sizeof(int),
                       STARPU_SCRATCH,   options->ws_worker,
                       STARPU_PRIORITY,    options->priority,
                       STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zgytlrb",
#endif
                       0);
}

/*   cl_zgytlr_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgytlrb_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int nAUV;
    MORSE_Complex64_t *AUV;
    MORSE_Complex64_t *AD = NULL;
    double *Ark;
  //  MORSE_Complex64_t *Dense;
    int lda;
    int ldu;
    int ldv;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;
    int maxrank;
    double tol;
    int compress_diag;
   int ntrian, nip, local_nt, NB; 
    AUV = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double*)STARPU_MATRIX_GET_PTR(descr[1]);
    MORSE_Complex64_t *work = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);

//    Dense = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    
    
    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &lda, &ldu, &ldv, &bigM, &m0, &n0, &seed, &maxrank, &tol, &compress_diag, &ntrian,  &nip, &local_nt,  &NB);

     int p=m0*local_nt+1;
     int q=n0*local_nt+1;
//    printf("\n nip:%d, ntrian:%d, NB:%d, local_nt:%d, p:%d, q:%d\n", nip, ntrian, NB, local_nt, p, q);
    AL4SAN_Matrix_Generation(&nip, &ntrian, work, &p, &q, &local_nt, &NB);
    
    MORSE_Complex64_t *AU = AUV;

   int nAU;
    nAU = nAUV/2;
    assert(ldu == ldv);
    size_t nelm_AU = (size_t)ldu * (size_t)nAU;
    MORSE_Complex64_t *AV = &(AUV[nelm_AU]);
    
    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,MORSE_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    HCORE_zgytlr( m, n,
                 AU,
                 AV,
                 AD,
                 Ark,
                 lda,
                 ldu,
                 ldv,
                 bigM, m0, n0, seed,
                 maxrank, tol,
                 compress_diag,
                 work
                 );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgytlrb, 3, cl_zgytlrb_cpu_func)
//#endif
