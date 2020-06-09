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

#include "runtime/starpu/runtime_codelets.h"
ZCODELETS_HEADER(cpy)
//extern int isPrecentageused;
//extern float precent2;

//CHAMELEON_CL_CB(zgytlr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N)
/*   MORSE_TASK_zgytlr - Generate a tile for random matrix. */

//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zcpy( const MORSE_option_t *options,
                       int m, int n,
                       const MORSE_desc_t *AUV,
                       const MORSE_desc_t *AUVnew,
                       const MORSE_desc_t *Ark,
                       int Am, int An,
                       int ldamUV,
                       int ldamUVnew,
                       int bigM, int m0, int n0, unsigned long long int seed,
                       int maxrank, MORSE_Complex64_t tol,
                       int compress_diag
                       )
{
    struct starpu_codelet *codelet = &cl_zcpy;
    //void (*callback)(void*) = options->profiling ? cl_zgytlr_callback : NULL;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;
    int nAUVnew = AUVnew->nb;
 
    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(AUV, Am, An);
    MORSE_ACCESS_W(AUVnew, Am, An);
    MORSE_ACCESS_W(Ark, Am, An);
    MORSE_END_ACCESS_DECLARATION;
    
    //printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);
    //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    //printf("\n nAUV:%d m:%d, n:%d nAUVnew:%d\n", nAUV, m, n, nAUVnew);

    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,    &m,                      sizeof(int),
                       STARPU_VALUE,    &n,                      sizeof(int),
                       STARPU_VALUE,    &nAUV,                      sizeof(int),
                       STARPU_VALUE,    &nAUVnew,                      sizeof(int),
                       STARPU_R,         RTBLKADDR(AUV, MORSE_Complex64_t, Am, An),
                       STARPU_W,         RTBLKADDR(AUVnew, MORSE_Complex64_t, Am, An), 
                       STARPU_R,         RTBLKADDR(Ark, double, Am, An),
                       STARPU_VALUE,  &ldamUV,                      sizeof(int),
                       STARPU_VALUE,  &ldamUVnew,                      sizeof(int),
                       STARPU_VALUE, &bigM,                      sizeof(int),
                       STARPU_VALUE,   &Am,                      sizeof(int),
                       STARPU_VALUE,   &An,                      sizeof(int),
                       STARPU_VALUE, &seed,   sizeof(unsigned long long int),
                       STARPU_VALUE,   &maxrank,                      sizeof(int),
                       STARPU_VALUE,   &tol,                      sizeof(MORSE_Complex64_t),
                       STARPU_VALUE,   &compress_diag,            sizeof(int),
                       STARPU_PRIORITY,    options->priority,
                       STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zcpy",
#endif
                       0);
}

/*   cl_zgytlr_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zcpy_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int nAUV;
    int nAUVnew;
    MORSE_Complex64_t *AUV;
    MORSE_Complex64_t *AUVnew;
    double *Ark;
    int ldamUV;
    int ldamUVnew;
    int bigM;
    int Am;
    int An;
    unsigned long long int seed;
    int maxrank;
    MORSE_Complex64_t tol;
    int compress_diag;
    
    AUV = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    AUVnew = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[1]);
    Ark = (double*)STARPU_MATRIX_GET_PTR(descr[2]);
    
    int _Ark;
    if(HICMA_get_always_fixed_rank() == 1){
        _Ark = HICMA_get_fixed_rank();
    } else {
        _Ark = Ark[0];
    }
    //printf("%s %s %d: _Ark:%d\n", __FILE__, __func__, __LINE__, _Ark); 
    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &nAUVnew, &ldamUV, &ldamUVnew, &bigM, &Am, &An, &seed, &maxrank, &tol, &compress_diag );
    
    
    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,MORSE_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
    HCORE_zcpy( m, n, nAUV, nAUVnew,
                 AUV,
                 AUVnew,
                 _Ark,
                 ldamUV,
                 ldamUVnew,
                 bigM, Am, An, seed,
                 maxrank, tol,
                 compress_diag
                 );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zcpy, 3, cl_zcpy_cpu_func)
//#endif
