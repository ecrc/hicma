




#include "morse.h"
#include "hicma.h"
#include "hicma_common.h"
#include "control/common.h"
#include "hicma_runtime_z.h"
#include "coreblas/lapacke.h"

#include "control/hicma_config.h"
#include <stdio.h>


#define A(p,q) A,  p, q  
#define GWASBLKLDD(A, k) A->get_blkldd( A, k )

/**
 *  gwas_psplgsy - adding regularizer.
 */
void compute_RHS_gene( int ntrian, int nip, int local_nt, int NB, MORSE_desc_t *A,
                    MORSE_sequence_t *sequence, MORSE_request_t *request )
{
    MORSE_context_t *morse;
    MORSE_option_t options;

    int p,q, counterp=0, counterq=0;
    int ldam;
    int tempmm, tempnn;

    morse = morse_context_self();
    if (sequence->status != MORSE_SUCCESS)
        return;

    //printf("\n ntrian:%d, nip:%d, local_nt:%d, NB:%d\n", ntrian, nip, local_nt, NB);
    RUNTIME_options_init(&options, morse, sequence, request);
  //for (q=0;q<A->nt;q=q+1){
       counterp=0;
        q=0;
    //printf("\n A->mt:%d\n",A->mt);
     for (p=0;p<A->mt;p=p+1){
       INSERT_TASK_RHS_gene(&options,
                        ntrian, nip, local_nt, NB, A(p,q), counterp, counterq);
      counterp=counterp+local_nt-1;
    }
       counterq=counterq+local_nt-1;
    // }
    RUNTIME_options_finalize(&options, morse);
}


