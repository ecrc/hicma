/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
 *                     All rights reserved.
 **/
/**
 * @file codelet_zgetrf.c
 *
 *  HiCMA codelets kernel
 *  HiCMA is a software package provided by King Abdullah University of Science and Technology (KAUST)
 *
 * @version 1.0.0
 * @author Kadir Akbudak
 * @date 2018-11-08
 **/

#include <sys/time.h>
#include <runtime/starpu/hicma_starpu.h>
#include <misc/auxdescutil.h>
#include <runtime/starpu/hicma_runtime_codelets.h>

ZCODELETS_HEADER(getrf_hcore)
//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zgetrf(const HICMA_option_t *options, int n, int nb,
                       const HICMA_desc_t *A, int Am, int An, int lda,
                       int iinfo)
{
    (void)nb;
    struct starpu_codelet *codelet = &cl_zgetrf_hcore;
    /*void (*callback)(void*) = options->profiling ? cl_zpotrf_callback : NULL;*/
    void (*callback)(void*) =  NULL;
    
    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_RW(A, Am, An);
    HICMA_END_ACCESS_DECLARATION;
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       //   STARPU_VALUE,    &uplo,                      sizeof(HICMA_enum),
                       STARPU_VALUE,    &n,                         sizeof(int),
                       STARPU_RW,        RTBLKADDR(A, HICMA_Complex64_t, Am, An),
                       STARPU_VALUE,    &lda,                       sizeof(int),
                       STARPU_VALUE,    &iinfo,                     sizeof(int),
                       STARPU_VALUE,    &Am,                        sizeof(int),
                       STARPU_VALUE,    &An,                        sizeof(int),
                       /* STARPU_SCRATCH,   options->ws_worker, */
                       STARPU_PRIORITY,  options->priority,
                       STARPU_CALLBACK,  callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zgetrf",
#endif
                       0);
}

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgetrf_hcore_cpu_func(void *descr[], void *cl_arg)
{
#ifdef HICMA_DISABLE_ALL_COMPUTATIONS
    return;
#endif
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    
    HICMA_enum uplo;
    int n;
    HICMA_Complex64_t *A;
    int lda;
    int iinfo;
    int info = 0;
    int Am;
    int An;
    
    
    A = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    
    starpu_codelet_unpack_args(cl_arg, &n, &lda, &iinfo, &Am, &An);
    if(HICMA_get_print_index() == 1){
        printf("%d+GETRF\t|AD(%d,%d)\n",HICMA_My_Mpi_Rank(), Am,An);
    }
    // if(HICMA_get_print_mat() == 1){
    //    printf("%d\tgetrf-input\n", __LINE__);
    //      _printmat(A, n, n, lda);
    // }
    //CORE_zpotrf(uplo, n, A, lda, &info);
    /* info = LAPACKE_dpotrf_work(
     LAPACK_COL_MAJOR,
     hicma_lapack_const(uplo),
     n, A, lda); */
    int *ipiv=(int*)malloc(hicma_max(1, hicma_min(n, n))*sizeof(int));
//printf("\n Noha1\n");    
//    _printmat_complex(A, n, n, lda);

    info = LAPACKE_zgetrf(LAPACK_COL_MAJOR, n, n,  A, lda, ipiv);
    
    free(ipiv);
//    _printmat_complex(A, n, n, lda);
 //   if(HICMA_get_print_mat() == 1){
           //printf("%d\tgetrf-output\n", __LINE__);
           //_printmat_complex(A, n, n, lda);
//      }
//printf("\n Noha2\n");
//    printf("%s %d EARLY EXIT\n", __FILE__, __LINE__);
//    fflush(stdout);
//    exit(-1);

    if(info != 0){
        fprintf(stderr, "%s\t|%d\t|Error in LAPACK getrf. info:%d. 1:n:%d 2:n:%d 4:lda:%d\n", __FILE__, __LINE__, info, n, n, lda);
        if(0) {
            int p,q;
            for(p=0; p<n;p++){
                for(q=0; q<10;q++){
                    printf("%.2e ", A[q*lda+p]);
                }
                printf("\n");
            }
            //getc(stdin);
        }
        exit(1);
    }
    if(HICMA_get_print_index() == 1 || HICMA_get_print_index_end() == 1){
        gettimeofday (&tvalAfter, NULL);
        printf("%d-GETRF\t|AD(%d,%d) N:%d LD:%d\t\t\t\t\tgetrf:%.4f\n",HICMA_My_Mpi_Rank(), Am, An,
               n, lda,
               (tvalAfter.tv_sec - tvalBefore.tv_sec)
               +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0
               );
    }
    if(0){
        char datebuf_start[128];
        time_t timer;
        struct tm* tm_info;
        struct tvalAfter;
        gettimeofday (&tvalAfter, NULL);
        time(&timer);
        tm_info = localtime(&timer);
        strftime(datebuf_start, 26, "%Y-%m-%d %H:%M:%S",tm_info);
        fprintf(stderr, "%s::%d\t", datebuf_start, Am);
    }
    
}

#ifdef CHAMELEON_USE_MAGMA
static void cl_zpotrf_hcore_cuda_func(void *descr[], void *cl_arg)
{
    cudaStream_t stream[2], currentt_stream;
    HICMA_enum uplo;
    int n;
    cuDoubleComplex *A;
    /* cuDoubleComplex *hA; */
    int lda;
    int iinfo;
    int info = 0;
    
    A  = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[0]);
    starpu_codelet_unpack_args(cl_arg, &uplo, &n, &lda, &iinfo);
    
    /* /\* */
    /*  *  hwork => nb*nb */
    /*  *\/ */
    /* hA = (cuDoubleComplex *)STARPU_MATRIX_GET_PTR(descr[1]); */
    
    /*      stream[0] = starpu_cuda_get_local_stream(); */
    /*      if ( cudaStreamCreate( stream+1 ) != CUDA_SUCCESS ){ */
    /*          fprintf(stderr, "Error while creating stream in codelet_zpotrf\n"); */
    /*          exit(-1); */
    /*      } */
    
    CUDA_zpotrf( uplo, n, A, lda, &info);
    
    cudaThreadSynchronize();
    /*      cudaStreamDestroy( stream[1] ); */
    
    return;
}
#endif
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
#if defined CHAMELEON_USE_MAGMA
CODELETS(zgetrf_hcore, 1, cl_zgetrf_hcore_cpu_func, cl_zgetrf_hcore_cuda_func, 0)
#else
CODELETS_CPU(zgetrf_hcore, 1, cl_zgetrf_hcore_cpu_func)
#endif
//#endif
