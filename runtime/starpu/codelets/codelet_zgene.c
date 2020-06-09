/**
 * @copyright (c) 2019 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zgene.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 0.1.1
 * @author Rabab Alomairy
 * @date 2019-09-08
 * @precisions normal z -> c d s
 **/

#include "morse.h"
#include "hicma.h"
#include "hicma_common.h"
#include "runtime/starpu/chameleon_starpu.h"
//#include "runtime/starpu/include/runtime_codelet_z.h"

#include <sys/time.h>

#include "runtime/starpu/runtime_codelets.h"

ZCODELETS_HEADER(gene)

void AL4SAN_Matrix_Generation(int *nip, int *ntrian, MORSE_Complex64_t *zz, int *q, int *p, int *local_nt, int *nb);


void INSERT_TASK_gene( const MORSE_option_t *options,
                        int ntrian, int nip, int local_nt, int NB, const MORSE_desc_t *A, int Am, int An, int counterp, int counterq)
{

    struct starpu_codelet *codelet = &cl_zgene;
    //void (*callback)(void*) = options->profiling ? cl_srls_callback : NULL;

    MORSE_BEGIN_ACCESS_DECLARATION;
    MORSE_ACCESS_W(A, Am, An);
    MORSE_END_ACCESS_DECLARATION;
//    printf("\n local_nt:%d, Am:%d, An:%d\n", local_nt, Am, An);
    starpu_insert_task(
        starpu_mpi_codelet(codelet),
        STARPU_VALUE,    &ntrian,                      sizeof(int),
        STARPU_VALUE,    &nip,                      sizeof(int),
        STARPU_W,         RTBLKADDR(A, CHAMELEON_Complex64_t, Am, An),
        STARPU_VALUE,  &Am,                      sizeof(int),
        STARPU_VALUE, &An,                      sizeof(int),
        STARPU_VALUE,   &local_nt,                      sizeof(int),
        STARPU_VALUE,   &NB,                      sizeof(int),
        STARPU_VALUE,   &counterp,                      sizeof(int),
        STARPU_VALUE,   &counterq,                      sizeof(int),
        STARPU_PRIORITY,    options->priority,
        //STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
        STARPU_NAME, "gene",
#endif
        0);
}


/*   cl_splgsy_cpu_func - Generate a tile for random symmetric (positive definite if 'bump' is large enough) matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgene_cpu_func(void *descr[], void *cl_arg)
{
    int ntrian;
    int nip;
    MORSE_Complex64_t *A;
    int p;
    int q;
    int local_nt;
    int NB;
    int counterp, counterq;
    A = (MORSE_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &ntrian, &nip, &p, &q, &local_nt, &NB, &counterp, &counterq);
//    p=p+1;
//    q=q+1;
//    p=p+counterp;
//    q=q+counterq;
     p=p*local_nt+1;
     q=q*local_nt+1;
    //printf("\n local_nt:%d, p:%d, q:%d\n", local_nt, p, q);
    AL4SAN_Matrix_Generation(&nip, &ntrian, A, &p, &q, &local_nt, &NB);


            //printf("%d\tgenetation-output\n", __LINE__);
           //_printmat_complex(A, NB, NB, NB);
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgene, 1, cl_zgene_cpu_func)
