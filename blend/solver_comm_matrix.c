/**
 *
 * @file solver_matrix_gen.c
 *
 * PaStiX solver structure generation function.
 *
 * @copyright 1998-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Tony Delarue
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#define _GNU_SOURCE 1
#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>

#include "common.h"
#include "pastix/order.h"
#include "cost.h"
#include "symbol.h"
#include "queue.h"
#include "solver.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_blend
 *
 * @brief Initialize the solver matrix structure
 *
 * This function takes all the global preprocessing steps: the symbol matrix,
 * and the resul of the simulation step to generate the solver matrix that hold
 * only local information to each PaStiX process.
 *
 *******************************************************************************
 *
 * @param[in] clustnum
 *          The index of the local PaStiX process.
 *
 * @param[inout] solvmtx
 *          On entry, the allocated pointer to a solver matrix structure.
 *          On exit, this structure holds alls the local information required to
 *          perform the numerical factorization.
 *
 * @param[in] symbmtx
 *          The global symbol matrix structure.
 *
 * @param[in] ordeptr
 *          The ordering structure.
 *
 * @param[in] simuctrl
 *          The information resulting from the simulation that will provide the
 *          data mapping, and the order of the task execution for the static
 *          scheduling.
 *
 * @param[in] ctrl
 *          The blend control structure that contains extra information
 *          computed during the analyze steps and the parameters of the analyze
 *          step.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS if success.
 * @retval PASTIX_ERR_OUTOFMEMORY if one of the malloc failed.
 *
 *******************************************************************************/
void
solverCommunicationMatrix( const char *dirtemp,
                           const SolverMatrix *solvmtx, int solv_seq )
{
    SolverCblk  *wcblk, *rcblk;
    SolverBlok  *blok;
    pastix_int_t i, k, n, m;
    pastix_int_t thrd_rid, thrd_wid, proc_rid, proc_wid, wsize;
    char *fname;
    FILE *f;

    if ( solv_seq ) {
        asprintf( &fname, "%s/cost_matrix.csv", dirtemp );
    }
    else {
        asprintf( &fname, "%s/cost_matrix_%ld.csv",
                  dirtemp, (long)(solvmtx->clustnum) );
    }
    f = fopen( fname, "w" );
    free(fname);
    fprintf( f, "# Proc Reader; Thread Reader; Proc Writer; Thread Writer; Cblk index; Nbr of elements written\n" );

    rcblk = solvmtx->cblktab;
    for (i=0; i<solvmtx->cblknbr; i++, rcblk++) {
        thrd_rid = rcblk->threadid;
        proc_rid = rcblk->ownerid;

        blok = rcblk->fblokptr;

        if ( rcblk->cblktype & CBLK_FANIN ) {
            /* Nothing to do. Done earlier */
            proc_rid = solvmtx->clustnum;
            thrd_wid = rcblk->threadid;
            proc_wid = rcblk->ownerid;

            assert( thrd_wid == -1 );
            assert( proc_rid != proc_wid );

            wsize = cblk_colnbr( rcblk ) * rcblk->stride;

            fprintf( f, "%ld;%ld;%ld;%ld;%ld;%ld\n",
                     (long)proc_rid, (long)-1,
                     (long)proc_wid, (long)-1,
                     (long)rcblk->gcblknum, wsize );
            continue;
        }

        if ( rcblk->cblktype & CBLK_RECV ) {
            /* We just have one receive */
            wcblk = solvmtx->cblktab + blok->fcblknm;

            thrd_wid = wcblk->threadid;
            proc_wid = wcblk->ownerid;

            assert( thrd_rid == -1 );
            assert( thrd_wid != -1 );
            assert( !(wcblk->cblktype & (CBLK_RECV | CBLK_FANIN)) );
            assert( rcblk->stride <= wcblk->stride );
            assert( cblk_colnbr( rcblk ) <= cblk_colnbr( wcblk ) );

            wsize = cblk_colnbr( rcblk ) * rcblk->stride;

            /* Already printed via FANIN */
            /* fprintf( f, "%ld;%ld;%ld;%ld;%ld;%ld\n", */
            /*          (long)proc_rid, (long)thrd_rid, */
            /*          (long)proc_wid, (long)thrd_wid, */
            /*          (long)wcblk->gcblknum, wsize ); */
            continue;
        }

        /* Skip diagonal block which does not generate contribution */
        blok++;
        for( ; blok < rcblk[1].fblokptr; blok++ ) {
            wcblk = solvmtx->cblktab + blok->fcblknm;

            thrd_wid = wcblk->threadid;
            proc_wid = wcblk->ownerid;

            /* Receive should never be on the write side */
            assert( !(wcblk->cblktype & CBLK_RECV) );

            k  = cblk_colnbr( rcblk );
            n  = blok_rownbr( blok );
            m  = rcblk->stride;
            m -= (rcblk->cblktype & CBLK_LAYOUT_2D) ? blok->coefind / k : blok->coefind;

            wsize = m * n;

            if ( wcblk->cblktype & CBLK_FANIN ) {
                assert( thrd_wid == -1 );
                assert( thrd_rid != -1 );
                assert( proc_rid != proc_wid );

                fprintf( f, "%ld;%ld;%ld;%ld;%ld;%ld\n",
                         (long)proc_rid, (long)thrd_rid,
                         (long)proc_rid, (long)-2,
                         (long)wcblk->gcblknum, wsize );
            }
            else {
                assert( thrd_rid != -1 );
                assert( thrd_wid != -1 );

                if ( (thrd_rid != thrd_wid) )
                {
                    assert( proc_rid == proc_wid );
                    fprintf( f, "%ld;%ld;%ld;%ld;%ld;%ld\n",
                             (long)proc_rid, (long)thrd_rid,
                             (long)proc_wid, (long)thrd_wid,
                             (long)wcblk->gcblknum, wsize );
                }
            }
        }
    }

    fclose(f);
}
