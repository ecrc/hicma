/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/

/**
 * @file codelet_zplrnt.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Ali Charara
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

/*
 * @copyright (c) 2009-2014 The University of Tennessee and The University
 *                          of Tennessee Research Foundation.
 *                          All rights reserved.
 * @copyright (c) 2012-2016 Inria. All rights reserved.
 * @copyright (c) 2012-2014 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordeaux. All rights reserved.
 */

/**
 *
 * @file codelet_zplrnt.c
 *
 * @copyright 2009-2014 The University of Tennessee and The University of
 *                      Tennessee Research Foundation. All rights reserved.
 * @copyright 2012-2016 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 ***
 *
 * @brief Chameleon zplrnt StarPU codelet
 *
 * @version 1.0.0
 * @comment This file has been automatically generated
 *          from Plasma 2.5.0 for MORSE 1.0.0
 * @author Piotr Luszczek
 * @author Pierre Lemarinier
 * @author Mathieu Faverge
 * @author Emmanuel Agullo
 * @author Cedric Castagnede
 * @date 2010-11-15
 * @precisions normal z -> c d s
 *
 */

#include <hicma.h>
#include <include/hicma_runtime_z.h>
#include <runtime/starpu/hicma_starpu.h>

ZCODELETS_HEADER(plrnt)

#define Rnd64_A 6364136223846793005ULL
#define Rnd64_C 1ULL
#define RndF_Mul 5.4210108624275222e-20f
#define RndD_Mul 5.4210108624275222e-20

#if defined(PRECISION_z) || defined(PRECISION_c)
#define NBELEM   2
#else
#define NBELEM   1
#endif

/*   HICMA_TASK_zplrnt - Generate a tile for random matrix. */

void HICMA_TASK_zplrnt(const HICMA_option_t *options,
                             int m, int n, const HICMA_desc_t *A, int Am, int An, int lda,
                             int bigM, int m0, int n0, unsigned long long int seed) {

    struct starpu_codelet *codelet = &cl_zplrnt;
    void (*callback)(void *) = NULL;
//    void (*callback)(void*) = options->profiling ? cl_zplrnt_callback : NULL;

    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(A, Am, An);
    HICMA_END_ACCESS_DECLARATION;

    starpu_insert_task(
            starpu_mpi_codelet(codelet),
            STARPU_VALUE, &m, sizeof(int),
            STARPU_VALUE, &n, sizeof(int),
            STARPU_W, RTBLKADDR(A, HICMA_Complex64_t, Am, An),
            STARPU_VALUE, &lda, sizeof(int),
            STARPU_VALUE, &bigM, sizeof(int),
            STARPU_VALUE, &m0, sizeof(int),
            STARPU_VALUE, &n0, sizeof(int),
            STARPU_VALUE, &seed, sizeof(unsigned long long int),
            STARPU_PRIORITY, options->priority,
            STARPU_CALLBACK, callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
            STARPU_NAME, "zplrnt",
#endif
            0);
}

static unsigned long long int Rnd64_jump(unsigned long long int n, unsigned long long int seed) {
    unsigned long long int a_k, c_k, ran;
    int i;

    a_k = Rnd64_A;
    c_k = Rnd64_C;

    ran = seed;
    for (i = 0; n; n >>= 1, ++i) {
        if (n & 1)
            ran = a_k * ran + c_k;
        c_k *= (a_k + 1);
        a_k *= a_k;
    }

    return ran;
}

/*   cl_zplrnt_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)

static void cl_zplrnt_cpu_func(void *descr[], void *cl_arg) {
    int m;
    int n;
    HICMA_Complex64_t *A;
    int lda;
    int bigM;
    int m0;
    int n0;
    unsigned long long int seed;

    A = (HICMA_Complex64_t *) STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &m, &n, &lda, &bigM, &m0, &n0, &seed);


    //  CORE_zplrnt - Generate a tile for random matrix.
    HICMA_Complex64_t *tmp = A;
    int64_t i, j;
    unsigned long long int ran, jump;

    jump = (unsigned long long int) m0 + (unsigned long long int) n0 * (unsigned long long int) bigM;

    for (j = 0; j < n; ++j) {
        ran = Rnd64_jump(NBELEM * jump, seed);
        for (i = 0; i < m; ++i) {
            *A = 0.5f - ran * RndF_Mul;
            ran = Rnd64_A * ran + Rnd64_C;
#if defined(PRECISION_z) || defined(PRECISION_c)
            *tmp += I*(0.5f - ran * RndF_Mul);
        ran   = Rnd64_A * ran + Rnd64_C;
#endif
            tmp++;
        }
        tmp += lda - i;
        jump += bigM;
    }
}




#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zplrnt, 1, cl_zplrnt_cpu_func)
