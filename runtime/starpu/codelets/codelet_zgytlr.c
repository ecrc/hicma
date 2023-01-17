/**
 * @copyright (c) 2017-2022 King Abdullah University of Science and Technology (KAUST).
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

#include <runtime/starpu/hicma_starpu.h>
#include <hcore_z.h>
#include <runtime/starpu/hicma_runtime_codelets.h>
#include <hicma.h>

ZCODELETS_HEADER(gytlr)

//HICMA_CHAM_CL_CB(zgytlr,        starpu_matrix_get_nx(task->handles[0]), starpu_matrix_get_ny(task->handles[0]), 0,                                                M*N)
/*   HICMA_TASK_zgytlr - Generate a tile for random matrix. */

//#if defined(HICMA_COMPLEX)
void HICMA_TASK_zgytlr( const HICMA_option_t *options,
                       int m, int n,
                       const HICMA_desc_t *AUV,
                       const HICMA_desc_t *Ark,
                       int Am, int An,
                       int lda,
                       int ldu,
                       int ldv,
                       int bigM, int m0, int n0, unsigned long long int seed,
                       int maxrank, double tol,
                       int compress_diag,
                       const HICMA_desc_t *Dense
                       )
{
    struct starpu_codelet *codelet = &cl_zgytlr;
    //void (*callback)(void*) = options->profiling ? cl_zgytlr_callback : NULL;
    void (*callback)(void*) = NULL;
    int nAUV = AUV->nb;
    
    HICMA_BEGIN_ACCESS_DECLARATION;
    HICMA_ACCESS_W(AUV, Am, An);
    HICMA_ACCESS_W(Ark, Am, An);
    HICMA_ACCESS_R(Dense, Am, An);
    HICMA_END_ACCESS_DECLARATION;
    
    //printf("%s:%d: Am:%d An:%d lda:%d bigM:%d m0:%d n0:%d\n ", __FILE__, __LINE__, Am, An, lda, bigM, m0, n0);
    //printf("%s %d: Am:%d An:%d ADm:%d ADn:%d ptr:%p\n", __func__, __LINE__, Am, An, ADm, ADn, ptr);
    //printf("\n nAUV:%d m:%d, n:%d\n", nAUV, m, n);
    starpu_insert_task(
                       starpu_mpi_codelet(codelet),
                       STARPU_VALUE,    &m,                      sizeof(int),
                       STARPU_VALUE,    &n,                      sizeof(int),
                       STARPU_VALUE,    &nAUV,                      sizeof(int),
                       STARPU_W,         RTBLKADDR(AUV, HICMA_Complex64_t, Am, An),
                       STARPU_W,         RTBLKADDR(Ark, double, Am, An),
                       STARPU_R,         RTBLKADDR(Dense, HICMA_Complex64_t, Am, An),
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
                       STARPU_PRIORITY,    options->priority,
                       STARPU_CALLBACK,    callback,
#if defined(CHAMELEON_CODELETS_HAVE_NAME)
                       STARPU_NAME, "hcore_zgytlr",
#endif
                       0);
}

#undef CBLAS_SADDR
#define CBLAS_SADDR( _val_ ) &(_val_)
extern int print_index;
extern int store_only_diagonal_tiles;
extern int global_check;
extern int print_mat;
int use_rsdd = 1; 
int zgytlr_print_index = 0;
int zgytlr_print_index_end = 0;
void HCORE_zgytlr( int m, int n, /*dimension of squareAD*/
        HICMA_Complex64_t *AU,
        HICMA_Complex64_t *AV,
        HICMA_Complex64_t *AD,
        double *Ark,
        int lda,
        int ldu,
        int ldv,
        int bigM, int ii, int jj, unsigned long long int seed,
        int maxrank, double tol, int compress_diag,
        HICMA_Complex64_t *Dense
        )
{
    int64_t i, j;
    //printf("m:%d n:%d bigM:%d m0:%d n0:%d\n", m, n, bigM, m0, n0);
    struct timeval tvalBefore, tvalAfter;  // removed comma
    gettimeofday (&tvalBefore, NULL);
    if(Ark[0] >= 1 && Ark[0]<maxrank) {
        maxrank = (int)Ark[0];
    }
    if(print_index || zgytlr_print_index){
        fprintf(stderr, "%d+GYTLR\t|(%d,%d) irk:%g mxrk:%d m:%d n:%d lda:%d ldu:%d ldv:%d\n",HICMA_My_Mpi_Rank(), ii, jj, Ark[0], maxrank, m, n, lda, ldu, ldv);
    }

    int shape[2];
    int rank = 0;

    STARSH_blrf* blrf =  HICMA_get_starsh_format();
    STARSH_cluster *RC = blrf->row_cluster, *CC = RC;
    void *RD = RC->data, *CD = RD;

    HICMA_Complex64_t *saveAD;
                                                    // allocate space for dense tile
    if((ii != jj && store_only_diagonal_tiles == 1) // if tile is off diagonal and
                                                    // and only diagonal tiles are stored in a tall and skinny matrix
                                                    // store_only_diagonal_tiles is here because
                                                    // code can use AD to store dense tiles
            ||                                      // OR
            compress_diag == 1) {                   // diagonals are also compressed so AD may not used perhaps
        saveAD = AD;
        //AD = malloc(sizeof(double) * m * n);
        AD = malloc(sizeof(HICMA_Complex64_t) * lda * n);
        assert(m==lda);
    }

    if(print_mat){
        printf("%d\tgytlr-UV-input\n", __LINE__);
        //_printmat_complex(AD, m, n, lda);
        //_printmat_complex(AU, m, rank, ldu);
        //_printmat_complex(AV, ldv, rank, ldv);
        _printmat_complex(Dense, m, n, lda);
    }

  blrf->problem->kernel(m, n, RC->pivot+RC->start[ii], CC->pivot+CC->start[jj],
            RD, CD, AD, lda);

    char chall = 'A';
    if(global_check ==1){ /** copy from Dense to AD */
        zlacpy_(&chall,
                &m, &n, AD, &lda, Dense, &lda);
    }
    if(ii==jj){
       Ark[0]==m;
      return;
      }
     
    //starsh_dense_zlrrsdd(m, n, AD, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);
    //lapack_int LAPACKE_zgesvd( int matrix_layout, char jobu, char jobvt, lapack_int m, lapack_int n, lapack_complex_double* a, lapack_int lda, double* s, lapack_complex_double* u, lapack_int ldu, lapack_complex_double* vt, lapack_int ldvt, double* superb );
    //https://software.intel.com/en-us/node/521150
    struct timeval tvalBefore_svd, tvalAfter_svd;  // removed comma
    gettimeofday (&tvalBefore_svd, NULL);
    int info = 0;
    if(use_rsdd) {
        HICMA_Complex64_t  *work;
        int *iwork;
        int oversample = 10;
        int mn = m;
        int mn2 = maxrank+oversample;
        if(mn2 > mn)
            mn2 = mn;
        // Get size of temporary arrays
        size_t lwork = n, lwork_sdd = (4*mn2+7)*mn2;
        if(lwork_sdd > lwork)
            lwork = lwork_sdd;
        lwork += (size_t)mn2*(2*n+m+mn2+1);
        size_t liwork = 8*mn2;
        // Allocate temporary arrays
        //STARSH_MALLOC(iwork, liwork);
        iwork = malloc(sizeof(*iwork) * liwork);
        if(iwork == NULL) {
            fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! liwork:%d", __FILE__, __func__, __LINE__, liwork);
            exit(-1);
        }
        //STARSH_MALLOC(work, lwork);
        work = malloc(sizeof(*work) * lwork);
        if(work == NULL) {
            fprintf(stderr, "%s %s %d:\t Allocation failed. No memory! lwork:%d", __FILE__, __func__, __LINE__, lwork);
            exit(-1);
        }
        //AD is m x n. AU and AV are m x maxrank and n x maxrank correspondingly
        //starsh_kernel_drsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        //starsh_dense_dlrrsdd(m, n, AD, AU, AV, &rank, maxrank, oversample, tol, work, lwork, iwork);
        double* svd_S = starsh_dense_zlrrsdd(m, n, AD, lda, AU, ldu,  AV, ldv, &rank, maxrank, oversample, tol, work, lwork, iwork);

        if(0){ /** Print singular values */
            char* str = calloc(10*n, sizeof(char));
            char* pstr = str;
            //pstr+=sprintf(pstr, "ii:%d jj:%d m:%d n:%d\n", ii, jj, m, n);
            int limsv=m;
            pstr+=sprintf(pstr, "%d,%d,%d:", ii, jj, limsv);
            for(int i = 0; i < limsv; i++){
                pstr+=sprintf(pstr, "%.2e ", svd_S[i]);
            }
            pstr+=sprintf(pstr, "\n");
            printf("%s", str);
            free(str);
        }

        free(work);
        free(iwork);

    } else {
        HICMA_Complex64_t* Dense_clone = calloc(m * lda, sizeof(HICMA_Complex64_t));
        zlacpy_(&chall,
                &m, &n, AD, &lda, Dense_clone, &lda);


        int nsv = m < n ? m : n;
        double* sigma = calloc(nsv, sizeof(double));
        double* svdsuperb = calloc(nsv, sizeof(double));
        HICMA_Complex64_t *__AV = calloc(nsv * ldv, sizeof(HICMA_Complex64_t));
        info = LAPACKE_zgesvd(LAPACK_COL_MAJOR, 'S', 'S', m, n, Dense_clone, lda, sigma, AU, ldu, __AV, ldv, svdsuperb);
        free(Dense_clone);
        //if(info == 0) the execution is successful.
        if(info < 0) { 
            printf("the %d-th parameter had an illegal value\n", info);
        }
        if(info > 0) {
            printf("if ?bdsqr did not converge, %d specifies how many superdiagonals of the intermediate bidiagonal form B did not converge to zero (see the description of the superb parameter for details).\n", info);
        }
        if(info != 0){
            fprintf(stderr,
                    "%s %d ERROR in LAPACKE_zgesvd() info=%d "
                    "1:m=%d, 2:n=%d, 3:Dense=%p, 4:lda=%d, 5:sigma=%p,"
                    "6:U=%p, 7:ldu=%d, 8:V=%p, 9:ldv:%d,"
                    "10:svdsuperb:%p"
                    "\n",
                    __FILE__, __LINE__, info,
                    m, n, Dense, lda, sigma,
                    AU, ldu, __AV, ldv,
                    svdsuperb);
            exit(-1);
        }

        for(i=0; i<nsv; i++){
            if(sigma[i] < tol){
                break;
            }
        }
        rank = i;
        if(print_mat){
            printf("Singular values: ");
            int i;
            for(i=0; i<nsv; i++){
                printf("%.2e ", sigma[i]);
            }
            printf("\n");
        }
        if(rank == nsv){ //means that tile is dense.
            rank = nsv;
            fprintf(stderr, "%s %s %d: Dense off-diagonal block (%d,%d)\n", __FILE__, __func__, __LINE__, ii, jj);
            exit(0);
        }
        int k;
        for(k = 0; k < rank; k++){
            double ddiagval = sigma[k];
            HICMA_Complex64_t zdiagval = ddiagval+0.0 * _Complex_I;
            /*printf("k:%d %+.4f%+.4fi\n", k, creal(zdiagval), cimag(zdiagval));*/
            cblas_zscal(n, CBLAS_SADDR(zdiagval), &__AV[k], ldv);
        } 
        LAPACKE_zge_trans(LAPACK_COL_MAJOR, rank, n,
                __AV, ldv, AV, ldv);

        free(sigma);
        free(svdsuperb);
        free(__AV);
    }
    gettimeofday (&tvalAfter_svd, NULL);
    double time_svd = 
                (tvalAfter_svd.tv_sec - tvalBefore_svd.tv_sec)
                +(tvalAfter_svd.tv_usec - tvalBefore_svd.tv_usec)/1000000.0;

        if(rank == 0) rank = 1;
        Ark[0] = rank;
        if(store_only_diagonal_tiles == 1) {
            assert(AD != saveAD);
            free(AD);
        }
    if(print_mat){
        printf("%d\tgytlr-UV-output\n", __LINE__);
        /*_printmat_complex(AD, m, n, lda);*/
        _printmat_complex(AU, m, rank, ldu);
        _printmat_complex(AV, n, rank, ldv);
        /*_printmat_complex(Dense, m, n, lda);*/
    }
    if(print_index || zgytlr_print_index_end){
        gettimeofday (&tvalAfter, NULL);
        fprintf(stderr, "%d-GYTLR\t|(%d,%d) rk:%g m:%d n:%d lda:%d ldu:%d ldv:%d\t\t\t\t\tGYTLR: %.4f svd:%.4f\n",HICMA_My_Mpi_Rank(),ii,jj,Ark[0],m, n, lda, ldu, ldv,
                (tvalAfter.tv_sec - tvalBefore.tv_sec)
                +(tvalAfter.tv_usec - tvalBefore.tv_usec)/1000000.0,
                time_svd
               );
    }
    return;
}
/*   cl_zgytlr_cpu_func - Generate a tile for random matrix. */

#if !defined(CHAMELEON_SIMULATION)
static void cl_zgytlr_cpu_func(void *descr[], void *cl_arg)
{
    int m;
    int n;
    int nAUV;
    HICMA_Complex64_t *AUV;
    HICMA_Complex64_t *AD = NULL;
    double *Ark;
    HICMA_Complex64_t *Dense;
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
    
    AUV = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[0]);
    Ark = (double*)STARPU_MATRIX_GET_PTR(descr[1]);
    Dense = (HICMA_Complex64_t *)STARPU_MATRIX_GET_PTR(descr[2]);
    
    
    starpu_codelet_unpack_args(cl_arg, &m, &n, &nAUV, &lda, &ldu, &ldv, &bigM, &m0, &n0, &seed, &maxrank, &tol, &compress_diag );
    
    HICMA_Complex64_t *AU = AUV;

   int nAU;
    nAU = nAUV/2;
    assert(ldu == ldv);
    size_t nelm_AU = (size_t)ldu * (size_t)nAU;
    HICMA_Complex64_t *AV = &(AUV[nelm_AU]);
    
    //printf("(%d,%d)%d %s %d %d\n", m0/m,n0/n,HICMA_My_Mpi_Rank(), __func__, __LINE__, AD == Dense);
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
                 Dense
                 );
}
#endif /* !defined(CHAMELEON_SIMULATION) */

/*
 * Codelet definition
 */
CODELETS_CPU(zgytlr, 3, cl_zgytlr_cpu_func)
//#endif
