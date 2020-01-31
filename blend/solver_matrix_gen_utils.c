/**
 *
 * @file solver_matrix_gen_utils.c
 *
 * PaStiX solver structure generation functions to factorize
 * solver_matric_gen.c .
 *
 * @copyright 1998-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Tony Delarue
 * @author Pascal Henon
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @author Mathieu Faverge
 * @date 2020-01-26
 *
 * @addtogroup blend_dev_solver
 * @{
 *
 **/
#include "common.h"
#include "symbol.h"
#include "solver.h"
#include "elimintree.h"
#include "cost.h"
#include "cand.h"
#include "pastix/order.h"
#include "extendVector.h"
#include "simu.h"
#include "solver_matrix_gen_utils.h"

static inline int
is_value_inside( pastix_int_t value, pastix_int_t *array,
                 pastix_int_t begin, pastix_int_t  end )
{
    pastix_int_t idx = (begin + end)/2;
    pastix_int_t nbegin, nend;

    if( value == array[idx] )
        return 1;

    nbegin = (value > array[idx]) ? idx : begin;
    nend   = (value < array[idx]) ? idx : end;

    if( (nbegin == begin) && (nend == end) )
        return 0;

    return is_value_inside(value, array, nbegin, nend );
}

static inline int
is_fblok_contributor( const symbol_blok_t *blok,
                      const symbol_blok_t *fblk,
                            pastix_int_t  *array,
                            pastix_int_t   size )
{
    if( is_value_inside(fblk->lcblknm, array, 0, size) )
    {
        return ( (fblk->frownum >= blok->frownum) &&
                 (fblk->lrownum <= blok->lrownum) );
    }
    return 0;
}
/**
 *******************************************************************************
 *
 * @brief Compute the bloknbr of a remote cblk
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix.
 *
 * @param[in] symbcblk
 *          The pointer to the current symbol cblk.
 *
 * @param[in] simuctrl
 *          The pointer to simuctrl structure.
 *
 * @param[in] ownerid
 *          Owner of the cblk.
 *
 * @return The number of bloks in the cblk
 *
 *******************************************************************************/
pastix_int_t
solvMatGen_remote_bloknbr( const symbol_matrix_t *symbmtx,
                           const symbol_cblk_t   *symbcblk,
                           const SimuCtrl        *simuctrl,
                                 int              ownerid )
{
    symbol_blok_t *symbblok, *fsymbblk;
    pastix_int_t   i, size, browidx, browmax;
    pastix_int_t   nbbloks, bloknum, bloknbr;
    pastix_int_t   brownbr = symbcblk[1].brownum - symbcblk[0].brownum;
    pastix_int_t  *ctrb_cblk;

    MALLOC_INTERN(ctrb_cblk, brownbr, pastix_int_t);
    /* Get the number of contributors with the same ownerid */
    size = 0;
    for ( i=0; i<brownbr; i++)
    {
        bloknum  = symbmtx->browtab[symbcblk->brownum + i ];
        symbblok = symbmtx->bloktab + bloknum;
        if( ownerid == simuctrl->bloktab[bloknum].ownerclust ) {
            ctrb_cblk[size] = symbblok->lcblknm;
            size++;
        }
    }
    assert(size <= brownbr);
    ctrb_cblk = (pastix_int_t*)realloc( ctrb_cblk, size * sizeof(pastix_int_t) );

    /* Compute the bloknbr of the cblk */
    nbbloks  = symbcblk[1].bloknum - symbcblk[0].bloknum - 1;
    bloknbr  = 1;

    /* The diagonal blok is already in */
    symbblok = symbmtx->bloktab + symbcblk->bloknum + 1;
    for( i=0; i<nbbloks; i++, symbblok++ )
    {
        browidx = symbmtx->cblktab[symbblok->fcblknm].brownum;
        browmax = symbmtx->cblktab[symbblok->fcblknm + 1].brownum;

        fsymbblk = symbmtx->bloktab + symbmtx->browtab[ browidx ];
        /* Check all regarding bloks */
        while ( (fsymbblk < symbblok) && (browidx < browmax) )
        {
            if( is_fblok_contributor(symbblok, fsymbblk, ctrb_cblk, size) &&
               (ownerid == simuctrl->bloktab[ fsymbblk - symbmtx->bloktab ].ownerclust) ) {
                bloknbr++;
                break;
            }

            browidx++;
            fsymbblk = symbmtx->bloktab + symbmtx->browtab[ browidx ];
        }
    }
    memFree_null(ctrb_cblk);
    return bloknbr;
}

/**
 *******************************************************************************
 *
 * @brief Compute the bloknbr of a remote cblk
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix.
 *
 * @param[in] symbcblk
 *          The pointer to the current symbol cblk.
 *
 * @param[in] simuctrl
 *          The pointer to simuctrl structure.
 *
 * @param[inout] solvcblk
 *          The pointer to current SolverCblk. The fcolnum and lcolnum will be
 *          updated.
 *
 *  @param[inout] solvblok
 *          The pointer to current solverBlok.
 *
 * @param[in] lcblknum
 *         Local cblk index in cblktab
 *
 * @param[in] layout2D
 *          Is layout2D ?
 *
 * @param[in] ownerid
 *          Owner of the cblk.
 *
 * @return The number of bloks in the cblk
 *
 *******************************************************************************/
pastix_int_t
solvMatGen_remote_bloks( const symbol_matrix_t *symbmtx,
                         const symbol_cblk_t   *symbcblk,
                         const SimuCtrl        *simuctrl,
                               SolverCblk      *solvcblk,
                               SolverBlok      *solvblok,
                               pastix_int_t     lcblknum,
                               pastix_int_t     layout2D,
                               int              ownerid )
{
    symbol_blok_t *symbblok, *fsymbblk;
    SolverBlok    *fblok = solvblok;
    pastix_int_t   i, size, browidx, browmax, brownbr;
    pastix_int_t   nbbloks, bloknum;
    pastix_int_t   fcolnum, lcolnum;
    pastix_int_t   frownum, lrownum;
    pastix_int_t   sfrownum, slrownum;
    pastix_int_t   nbcols, stride = 0;
    pastix_int_t  *ctrb_cblk;
    int            exist;

    symbblok = symbmtx->bloktab + symbcblk->bloknum;
    brownbr  = symbcblk[1].brownum - symbcblk[0].brownum;
    browidx  = symbmtx->cblktab[symbblok->fcblknm].brownum;
    fcolnum  = symbmtx->cblktab[symbmtx->cblknbr - 1].lcolnum;
    lcolnum  = 0;

    /* First loop : get accurate fcolnum and lcolnum */
    MALLOC_INTERN(ctrb_cblk, brownbr, pastix_int_t);
    size = 0;
    for( i=0; i<brownbr; i++, browidx++ ) {
        fsymbblk = symbmtx->bloktab + symbmtx->browtab[ browidx ];
        bloknum  = fsymbblk - symbmtx->bloktab;
        if( ownerid != simuctrl->bloktab[bloknum].ownerclust ) {
            continue;
        }

        ctrb_cblk[size] = fsymbblk->lcblknm;
        size++;
        solvMatGen_get_rownum( symbmtx, fsymbblk, &frownum, &lrownum );
        fcolnum = pastix_imin( fcolnum, frownum );
        lcolnum = pastix_imax( lcolnum, lrownum );
    }
    assert(size <= brownbr);
    ctrb_cblk = (pastix_int_t*)realloc( ctrb_cblk, size * sizeof(pastix_int_t) );

    solvcblk->fcolnum = fcolnum;
    solvcblk->lcolnum = lcolnum;
    nbcols  = lcolnum - fcolnum + 1;
    nbbloks = symbcblk[1].bloknum - symbcblk[0].bloknum;
    for( i=0; i<nbbloks; i++, symbblok++ )
    {
        browidx  = symbmtx->cblktab[symbblok->fcblknm].brownum;
        browmax  = symbmtx->cblktab[symbblok->fcblknm + 1].brownum;

        fsymbblk = symbmtx->bloktab + symbmtx->browtab[ browidx ];
        sfrownum = symbblok->lrownum;
        slrownum = symbblok->frownum;
        exist    = 0;
        while( (fsymbblk < symbblok) && (browidx < browmax) )
        {
            if( is_fblok_contributor(symbblok, fsymbblk, ctrb_cblk, size) &&
               (ownerid == simuctrl->bloktab[fsymbblk - symbmtx->bloktab].ownerclust) ) {

                solvMatGen_get_rownum( symbmtx, fsymbblk, &frownum, &lrownum );
                sfrownum = pastix_imin( sfrownum, frownum );
                slrownum = pastix_imax( slrownum, lrownum );
                exist = 1;
            }
            browidx++;
            fsymbblk = symbmtx->bloktab + symbmtx->browtab[ browidx ];
        }
        if( exist ) {
            assert(sfrownum >= symbblok->frownum );
            assert(slrownum <= symbblok->lrownum );

            solvMatGen_init_blok( solvblok,
                                  lcblknum, -1,
                                  sfrownum, slrownum,
                                  stride, nbcols, layout2D );
            stride += (slrownum - sfrownum + 1);
            solvblok->gbloknm = -1;
            solvblok++;
        }
    }
    assert( (solvblok - fblok) ==  solvMatGen_remote_bloknbr(symbmtx, symbcblk, simuctrl, ownerid) );
    memFree_null(ctrb_cblk);

    return stride;
}


/**
 *******************************************************************************
 *
 * @brief Fill the local numbering arrays to compress the symbol information
 *        into solver.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix structure.
 *
 * @param[in] simuctrl
 *          The pointer to the simuctrl structure.
 *
 *  @param[inout] solvmtx
 *          Pointer to the solver matrix.
 *
 * @param[inout] cblklocalnum
 *          Local cblk infos.
 *
 * @param[inout] bloklocalnum
 *          Local blok infos.
 *
 * @param[inout] tasklocalnum
 *          Local tasks infos.
 *
 * @param[inout] fcbklocalnum
 *          Array to compute CBLK_FANIN infos.
 *
 * @param[inout] pcbklocalnum
 *          Array to compute CBLK_RECV infos.
 *
 *******************************************************************************/
void
solvMatGen_fill_localnums( const symbol_matrix_t *symbmtx,
                           const SimuCtrl  *simuctrl,
                           SolverMatrix    *solvmtx,
                           pastix_int_t    *cblklocalnum,
                           pastix_int_t    *bloklocalnum,
                           pastix_int_t    *tasklocalnum,
                           pastix_int_t    *fcbklocalnum,
                           pastix_int_t    *pcbklocalnum,
                           int            **recv_sources_ptr )
{
    pastix_int_t  *localindex;
    pastix_int_t  *countcluster;
    symbol_cblk_t *symbcblk;
    pastix_int_t   cblknum, fcbknum, brownum, brownbr;
    pastix_int_t   faninnbr, recvnbr;
    pastix_int_t   i, j, k, c, fc;
    pastix_int_t   flaglocal;
    pastix_int_t   clustnum = solvmtx->clustnum;
    int           *recv_sources = *recv_sources_ptr;

    /* Initialize the set of cluster candidates for each cblk */
    MALLOC_INTERN( countcluster, solvmtx->clustnbr, pastix_int_t );
    MALLOC_INTERN( localindex,   solvmtx->clustnbr, pastix_int_t );

    memset( localindex, 0, solvmtx->clustnbr * sizeof(pastix_int_t) );

    /*
     * Compute local number of tasks on each cluster
     */
    for ( i = 0; i < simuctrl->tasknbr; i++ ) {
        c = simuctrl->bloktab[ simuctrl->tasktab[i].bloknum ].ownerclust;

        tasklocalnum[i] = localindex[c];
        localindex[c]++;
    }
    solvmtx->tasknbr = localindex[clustnum];

    /*
     * Compute the local numbering of the fan-in and recv cblks on each cluster
     */

    /* Set the arrays to compute local informations */
    memset( localindex,   0,    solvmtx->clustnbr * sizeof( pastix_int_t ) );
    memset( pcbklocalnum, 0,    symbmtx->cblknbr  * sizeof( pastix_int_t ) );
    memset( fcbklocalnum, 0xff, symbmtx->cblknbr  * sizeof( pastix_int_t ) );

    cblknum  = 0;
    fcbknum  = 0;
    brownum  = 0;
    recvnbr  = 0;
    faninnbr = 0;
    symbcblk = symbmtx->cblktab;
    for ( i = 0; i < symbmtx->cblknbr; i++, symbcblk++ ) {
        brownbr = symbcblk[1].brownum - symbcblk[0].brownum;
        memset( countcluster, 0, solvmtx->clustnbr * sizeof( pastix_int_t ) );

        /*
         * The cblk is considered local if data are local, or if we store a
         * compressed copy for fanin
         */
        flaglocal = ( fcbklocalnum[i] != -1 ) || ( simuctrl->cblktab[i].owned );
        if ( !flaglocal ) {
            cblklocalnum[i] = -i - 1;
            continue;
        }

        /*
         * The cblk is local and we may receive remote information, let's:
         *    - compute the size of the compressed browtab
         *    - compute the number of remote fanin to be received for the update
         *    - compute the set of remote nodes sending a fanin
         *
         * Work on the incoming edges.
         */
        if ( ( simuctrl->cblktab[i].owned ) &&
             ( symbcblk[1].brownum > symbcblk[0].brownum ) )
        {
            for ( j = symbcblk[0].brownum; j < symbcblk[1].brownum; j++ ) {
                k = symbmtx->browtab[j];
                c = simuctrl->bloktab[k].ownerclust;

                assert( i == symbmtx->bloktab[k].fcblknm );

                /* Compute the amount of remote contributions */
                if ( c != clustnum ) {
                    countcluster[c]++;
                    brownbr--;
                }
            }
            assert( brownbr >= 0 );

            /* Compute the amount of different contributers in a cblk */
            for ( j=0; j<solvmtx->clustnbr; j++ )
            {
                /* Each time countcluster[j] != 0, a new sender can be registered*/
                if( countcluster[j] != 0 ){
                    pcbklocalnum[i]++;                              /* Amount of receptions for cblk i */
                    recv_sources[recvnbr] = j;                      /* Register source                 */
                    /* Duplicate the amount of for this owner */
                    localindex[clustnum] += solvMatGen_remote_bloknbr( symbmtx,symbcblk, simuctrl, j );
                    brownbr++;                                      /* One more blok will be in the browtab */
                    cblknum++;                                      /* Add one cblk                    */
                    recvnbr++;                                      /* Add one reception count         */
                }
            }
        }

        /*
         * If the cblk is local (fanin or classic), we add it to the list
         *
         * Work on the outgoing edges.
         */
        for ( j = symbcblk[0].bloknum; j < symbcblk[1].bloknum; j++ ) {
            c               = simuctrl->bloktab[j].ownerclust;
            bloklocalnum[j] = localindex[clustnum];
            localindex[clustnum]++;

            if ( c == clustnum ) {
                fcbknum = symbmtx->bloktab[j].fcblknm;
                fc      = simuctrl->bloktab[symbmtx->cblktab[fcbknum].bloknum].ownerclust;

                /* If the facing cblk isn't local, we need to have a local copy of it */
                if ( ( fcbklocalnum[fcbknum] == -1 ) && ( fc != c ) ) {
                    fcbklocalnum[fcbknum] = fc;
                    faninnbr++;
                }
            }
        }
        /*
         * If this is a fanin, we compute the size of the local browtab
         */
        if ( fcbklocalnum[i] != -1 ) {
            localindex[clustnum] -= ( symbcblk[1].bloknum - symbcblk[0].bloknum
                                    - solvMatGen_remote_bloknbr( symbmtx,symbcblk, simuctrl, clustnum ) );
            for ( j = symbcblk[0].brownum; j < symbcblk[1].brownum; j++ ) {
                k = symbmtx->browtab[j];
                c = simuctrl->bloktab[k].ownerclust;
                if ( c != clustnum ) {
                    brownbr--;
                }
            }
        }

        /* Store index of the current cblk */
        cblklocalnum[i] = cblknum;
        cblknum++;

        /* Update the brownum index */
        brownum += brownbr;
        assert( brownum <= symbcblk[1].brownum );
    }

    solvmtx->cblknbr = cblknum;
    solvmtx->bloknbr = localindex[clustnum];
    solvmtx->brownbr = brownum;

    /* Reallocate recv_sources tab to diminish it's size */
    if( recvnbr > 0 ){
        int *sources = (int *)realloc( recv_sources, recvnbr * sizeof(int) );
        *recv_sources_ptr = sources;
    }
    solvmtx->recvnbr  = recvnbr;
    solvmtx->faninnbr = faninnbr;

    memFree_null( localindex );
    memFree_null( countcluster );
}

/**
 *******************************************************************************
 *
 * @brief Reorder the browtab from the symbol structure in a distributed way.
 *        First stock the 1D blocks and then the 2D blocks.
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix structure.
 *
 * @param[in] symbcblk
 *          The pointer to the current symbol cblk.
 *
 *  @param[inout] solvmtx
 *          Pointer to the solver matrix.
 *
 * @param[inout] solvcblk
 *          The pointer to the current solver cblk.
 *
 * @param[inout] browtmp
 *          An array that we will fill with the local browtab infos.
 *
 * @param[in] cblklocalnum
 *          Local cblk infos.
 *
 * @param[in] bloklocalnum
 *          Local blok infos.
 *
 * @param[in] brownum
 *         Current brownum.
 *
 *******************************************************************************/
pastix_int_t
solvMatGen_reorder_browtab( const symbol_matrix_t *symbmtx,
                            symbol_cblk_t         *symbcblk,
                            SolverMatrix          *solvmtx,
                            SolverCblk            *solvcblk,
                            pastix_int_t          *browtmp,
                            pastix_int_t          *cblklocalnum,
                            pastix_int_t          *bloklocalnum,
                            pastix_int_t           brownum )
{
    pastix_int_t   brownbr;
    symbol_blok_t *symbblok;
    SolverBlok    *solvblok;
    SolverCblk    *browcblk;
    pastix_int_t   lcblknm, lbloknm;
    pastix_int_t   j2d, j1d, j, jmax;
    pastix_int_t   *b;

    brownbr = symbcblk[1].brownum - symbcblk[0].brownum;
    solvcblk->brown2d = solvcblk->brownum + brownbr;

    /* Nothing to do here */
    if ( !brownbr ) {
        return 0;
    }

    assert( brownbr <= symbmtx->browmax );
    memcpy( browtmp,
            symbmtx->browtab + symbcblk->brownum,
            brownbr * sizeof(pastix_int_t) );

    /*
     * j   is the index in the local browtab (~ postition of b)
     * j1d is the number of discovered 1d block in the browtab
     * j2d if the index of the first 2d block in the original tab
     * jmax is equal to brownbr
     * brownbr is updated to store the real number of brow (minus fanin)
     *
     * b is a pointer to the temporary copy of the subsection of the browtab
     * At the end of the first pass, if b[i] is negative, it has already been treated
     * with if equal to:
     *       -1, the block was a 1D block
     *       -2, the block belonged to a remote cblk
     *       -3, the block belonged to a local fanin (should not happen)
     * It is is positive, it's a 2D block that need to be pushed to the end of
     * the browtab in the second pass.
     */
    b = browtmp;
    j2d = -1;
    jmax = brownbr;

    /* First pass to copy 1D updates */
    for ( j=0, j1d=0; j < jmax; j++, b++ ) {
        /* Get the contributing block in the symbol */
        symbblok = symbmtx->bloktab + (*b);

        lcblknm = ( cblklocalnum == NULL ) ? symbblok->lcblknm : cblklocalnum[ symbblok->lcblknm ];

        /* If distant blok */
        if ( lcblknm < 0 ) {
            *b = -2;
            brownbr--;
            continue;
        }

        /* Get the local cblk which owns the block */
        browcblk = solvmtx->cblktab + lcblknm;

        /* Recv should never appear through cblklocalnum */
        assert( !(browcblk->cblktype & CBLK_RECV) );

        /* Fanin should not contribute to local data */
        if( browcblk->cblktype & CBLK_FANIN ) {
            *b = -3;
            brownbr--;
            continue;
        }

        /* Store the first non 1D index to not rediscover the begining, and skip 2d for now */
        if ( browcblk->cblktype & CBLK_TASKS_2D ) {
            j2d = ( j2d == -1 ) ? j : j2d;
            continue;
        }

        /* Find the SolvBlok corresponding to the SymbBlok */
        lbloknm = ( bloklocalnum == NULL ) ? *b : bloklocalnum[ *b ];
        solvblok = solvmtx->bloktab + lbloknm;

        assert( solvblok->lcblknm == lcblknm );
        assert( ( symbblok->frownum == solvblok->frownum ) &&
                ( symbblok->lrownum == solvblok->lrownum ) );

        solvmtx->browtab[brownum + j1d] = lbloknm;
        solvblok->browind = brownum + j1d;
        *b = -1;
        j1d++;
    }

    /* Store the index of the first 2D contribution in the array */
    assert( j1d <= brownbr );
    solvcblk->brown2d = solvcblk->brownum + j1d;

    /* Second pass to copy 2D updates to the end */
    if ( j2d != -1 ) {
        b = browtmp + j2d;
        for ( j = j2d; j < jmax; j++, b++ ) {
            symbblok = symbmtx->bloktab + ( *b );

            if ( *b < 0 ) {
                continue;
            }
            lcblknm = ( cblklocalnum == NULL ) ? symbblok->lcblknm : cblklocalnum[ symbblok->lcblknm ];
            assert( lcblknm >= 0 );

            /* Get the local cblk which owns the block */
            browcblk = solvmtx->cblktab + lcblknm;
            assert( (cblklocalnum == NULL) ||
                    (browcblk->ownerid == solvmtx->clustnum) );

            /* Find the SolvBlok corresponding to the SymbBlok */
            lbloknm = ( bloklocalnum == NULL ) ? *b : bloklocalnum[ *b ];
            solvblok = solvmtx->bloktab + lbloknm;

            assert( solvblok->lcblknm == lcblknm );
            assert( ( symbblok->frownum == solvblok->frownum ) &&
                    ( symbblok->lrownum == solvblok->lrownum ) );

            solvmtx->browtab[brownum + j1d] = lbloknm;
            solvblok->browind = brownum + j1d;
            j1d++;
        }
    }
    assert( j1d == brownbr );

    return brownbr;
}

struct args_ttsktab
{
    SolverMatrix   *solvmtx;
    const SimuCtrl *simuctrl;
    pastix_int_t   *tasklocalnum;
    pastix_int_t    clustnum;
};

/**
 *******************************************************************************
 *
 * @brief Fill the ttsktab for it's own thread.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          the context of the current thread
 *
 * @param[inout] args
 *          The parameter as specified in solverMatricGen*.
 *
 *******************************************************************************/
void
solvMatGen_fill_ttsktab( isched_thread_t *ctx, void *args )
{
    struct args_ttsktab *arg     = (struct args_ttsktab*)args;
    SolverMatrix   *solvmtx      = arg->solvmtx;
    const SimuCtrl *simuctrl     = arg->simuctrl;
    pastix_int_t   *tasklocalnum = arg->tasklocalnum;
    pastix_int_t    clustnum     = arg->clustnum;
    int             rank         = ctx->rank;
    SimuProc       *simuproc     = simuctrl->proctab
                               + ( simuctrl->clustab[clustnum].fprocnum + rank );
    pastix_int_t i;
    pastix_int_t priomin = PASTIX_INT_MAX;
    pastix_int_t priomax = 0;
    pastix_int_t ttsknbr = extendint_Size( simuproc->tasktab );
    pastix_int_t j, jloc, cblknum;

    solvmtx->ttsknbr[rank] = ttsknbr;
    if(ttsknbr > 0) {
        MALLOC_INTERN(solvmtx->ttsktab[rank], ttsknbr, pastix_int_t);
    }
    else {
        solvmtx->ttsktab[rank] = NULL;
    }

    for(i=0; i<ttsknbr; i++)
    {
        j = extendint_Read(simuproc->tasktab, i);
        if( tasklocalnum != NULL ){
            jloc = tasklocalnum[j];
        }
        else {
            jloc = j;
        }
        /* Only local cblks should appear in the tasktab */
        cblknum = solvmtx->tasktab[jloc].cblknum;
        assert( !(solvmtx->cblktab[ cblknum ].cblktype & CBLK_FANIN) );
        solvmtx->ttsktab[rank][i] = jloc;

        solvmtx->cblktab[cblknum].threadid = rank;

        priomax = pastix_imax( solvmtx->tasktab[jloc].prionum, priomax );
        priomin = pastix_imin( solvmtx->tasktab[jloc].prionum, priomin );
    }

#if defined(PASTIX_DYNSCHED)
    solvmtx->btree->nodetab[rank].priomin = priomin;
    solvmtx->btree->nodetab[rank].priomax = priomax;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Fill in ttsktab for it's own thread. Only for debugging factorization.
 *
 *******************************************************************************
 *
 * @param[in] ctx
 *          the context of the current thread
 *
 * @param[inout] args
 *          The parameter as specified in solverMatricGen*.
 *
 *******************************************************************************/
void
solvMatGen_fill_ttsktab_dbg( isched_thread_t *ctx, void *args )
{
    struct args_ttsktab *arg = (struct args_ttsktab*)args;

    pastix_int_t  i, j, jloc, size;
    SolverMatrix *solvmtx = arg->solvmtx;
    int           rank    = ctx->rank;
    int           nthread = ctx->global_ctx->world_size;
    pastix_int_t  tasknbr = solvmtx->tasknbr / nthread;
    pastix_int_t  priomin = PASTIX_INT_MAX;
    pastix_int_t  priomax = 0;

    size = (rank == nthread-1) ? (solvmtx->tasknbr - (nthread-1) * tasknbr) : tasknbr;
    solvmtx->ttsknbr[rank] = size;

    if(size > 0) {
        MALLOC_INTERN(solvmtx->ttsktab[rank], size, pastix_int_t);
    }
    else {
        solvmtx->ttsktab[rank] = NULL;
    }

    j = ((solvmtx->tasknbr - (nthread-1) * tasknbr) * rank);
    for(i=0; i < size; i++)
    {
        solvmtx->ttsktab[rank][i] = j;

        jloc = solvmtx->tasktab[j].cblknum;
        solvmtx->cblktab[jloc].threadid = rank;

        priomax = pastix_imax( solvmtx->tasktab[j].prionum, priomax );
        priomin = pastix_imin( solvmtx->tasktab[j].prionum, priomin );
        j++;
    }

#if defined(PASTIX_DYNSCHED)
    solvmtx->btree->nodetab[rank].priomin = priomin;
    solvmtx->btree->nodetab[rank].priomax = priomax;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Fill the tasktab.
 *
 *******************************************************************************
 *
 * @param[inout] solvmtx
 *          Pointer to the solver matrix.
 *
 * @param[in] isched
 *          The internal context to run multi-threaded functions.
 *
 * @param[in] simuctrl
 *          The pointer to the simuctrl structure.
 *
 * @param[in] tasklocalnum
 *          Local tasks infos.
 *
 * @param[in] cblklocalnum
 *          Local cblk infos.
 *
 * @param[in] bloklocalnum
 *          Local blok infos.
 *
 * @param[in] clustnum
 *          Rank of the MPI instance.
 *
 * @param[in] is_dbg
 *          Enable/disable the ttsktab debug generation.
 *
 *******************************************************************************/
void
solvMatGen_fill_tasktab( SolverMatrix   *solvmtx,
                         isched_t       *isched,
                         const SimuCtrl *simuctrl,
                         pastix_int_t   *tasklocalnum,
                         pastix_int_t   *cblklocalnum,
                         pastix_int_t   *bloklocalnum,
                         pastix_int_t    clustnum,
                         int             is_dbg )
{
    SimuTask    *simutask = simuctrl->tasktab;
    Task        *solvtask;
    pastix_int_t nbftmax  = 0;
    pastix_int_t tasknum  = 0;
    pastix_int_t i;

    MALLOC_INTERN( solvmtx->tasktab, solvmtx->tasknbr+1, Task );
    solvtask = solvmtx->tasktab;

    /* No local indices, this is a global solver */
    if ( tasklocalnum == NULL )
    {
        for(i=0; i<simuctrl->tasknbr; i++, simutask++)
        {
            nbftmax = pastix_imax( nbftmax, simutask->ftgtcnt );
            assert( tasknum == i );

            solvtask->taskid  = COMP_1D;
            solvtask->prionum = simutask->prionum;
            solvtask->cblknum = simutask->cblknum;
            solvtask->bloknum = simutask->bloknum;
            solvtask->ctrbcnt = simutask->ctrbcnt;

            tasknum++; solvtask++;
        }
    }
    else
    {
        for(i=0; i<simuctrl->tasknbr; i++, simutask++)
        {
            if( simuctrl->bloktab[ simutask->bloknum ].ownerclust == clustnum )
            {
                assert( tasknum == tasklocalnum[i] );

                solvtask->taskid  = COMP_1D;
                solvtask->prionum = simutask->prionum;
                solvtask->cblknum = cblklocalnum[ simutask->cblknum ];
                solvtask->bloknum = bloklocalnum[ simutask->bloknum ];
                solvtask->ctrbcnt = simutask->ctrbcnt;

                tasknum++; solvtask++;
            }
        }
    }
    assert(tasknum == solvmtx->tasknbr);

    /* One more task to avoid side effect */
    solvtask->taskid  = -1;
    solvtask->prionum = -1;
    solvtask->cblknum = solvmtx->cblknbr+1;
    solvtask->bloknum = solvmtx->bloknbr+1;
    solvtask->ctrbcnt = 0;

    solvmtx->nbftmax = nbftmax;

    /* Fill in the ttsktab arrays (one per thread) */
    MALLOC_INTERN(solvmtx->ttsknbr, solvmtx->bublnbr, pastix_int_t  );
    MALLOC_INTERN(solvmtx->ttsktab, solvmtx->bublnbr, pastix_int_t* );

    if( is_dbg ) {
        struct args_ttsktab args = { solvmtx, NULL, tasklocalnum, clustnum };
        isched_parallel_call( isched, solvMatGen_fill_ttsktab_dbg, &args );
    }
    else {
        struct args_ttsktab args = { solvmtx, simuctrl, tasklocalnum, clustnum };
        isched_parallel_call( isched, solvMatGen_fill_ttsktab, &args );
    }
}

/**
 *******************************************************************************
 *
 * @brief Compute the maximum area of the temporary buffer
 *        used during computation
 *
 * During this loop, we compute the maximum area that will be used as
 * temporary buffers, and statistics:
 *    - diagmax: Only for hetrf/sytrf factorization, this the maximum size
 *               of a panel of MAXSIZEOFBLOCKS width in a diagonal block
 *    - gemmmax: For all, this is the maximum area used to compute the
 *               compacted gemm on a CPU.
 *
 * Rk: This loop is not merged with the main block loop, since strides have
 * to be peviously computed.
 *
 *******************************************************************************
 *
 *  @param[inout] solvmtx
 *           Pointer to the solver matrix.
 *
 *******************************************************************************/
void
solvMatGen_max_buffers( SolverMatrix *solvmtx )
{
    SolverCblk  *solvcblk = solvmtx->cblktab;
    SolverBlok  *solvblok = solvmtx->bloktab;
    pastix_int_t gemmmax = 0;
    pastix_int_t offdmax = 0;
    pastix_int_t blokmax = 0;
    pastix_int_t gemmarea, offdarea, cblk_m, acc_m, i;

    for(i=0; i<solvmtx->cblknbr; i++, solvcblk++)
    {
        SolverBlok *lblok = solvcblk[1].fblokptr;
        pastix_int_t m = solvcblk->stride;
        pastix_int_t n = cblk_colnbr( solvcblk );
        pastix_int_t k = blok_rownbr( solvblok );

        m -= n;

        /*
         * Compute the surface of the off-diagonal block in a panel for
         * LDL^[th] factorizations
         */
        offdarea = m * n;
        offdmax = pastix_imax( offdmax, offdarea );

        /*
         * Compute the maximum area for 1d temporary workspace in GEMM
         */
        solvblok++;
        cblk_m = -1;
        acc_m  = 0;
        for( ; solvblok<lblok; solvblok++ ) {
            k = blok_rownbr( solvblok );

            /*
             * Temporary workspace for GEMM
             * m+1 to store the diagonal in case of GEMDM
             */
            if ( !(solvcblk->cblktype & CBLK_LAYOUT_2D) ) {
                gemmarea = (m+1) * k;
                gemmmax = pastix_imax( gemmmax, gemmarea );
            }

            /*
             * Max size for off-diagonal blocks for 2-terms version of the
             * 2D LDL
             */
            if ( solvcblk->cblktype & (CBLK_TASKS_2D | CBLK_COMPRESSED) ) {
                if ( solvblok->fcblknm == cblk_m ) {
                    acc_m += k;
                }
                else {
                    cblk_m = solvblok->fcblknm;
                    acc_m = k;
                }
                blokmax = pastix_imax( n * acc_m, blokmax );
            }
            m -= k;
        }
    }

    solvmtx->offdmax = offdmax;
    solvmtx->gemmmax = gemmmax;
    solvmtx->blokmax = blokmax;
}

/**
 *******************************************************************************
 *
 * @brief Mark blocks if they belong to the last supernode, or if they are
 * facing it for statistical purpose only.
 *
 * TODO : Should be improved by using the brow array in order to cover the
 *        blocks in front of the last cblk
 *
 *******************************************************************************
 *
 *  @param[inout] solvmtx
 *           Pointer to the solver matrix.
 *
 *******************************************************************************/
void
solvMatGen_stats_last( SolverMatrix *solvmtx )
{
#if defined(PASTIX_SUPERNODE_STATS)
    pastix_int_t i;
    SolverBlok  *solvblok = solvmtx->bloktab;

    for(i=0; i<solvmtx->bloknbr; i++, solvblok++ ) {
        SolverCblk *fcblk = solvmtx->cblktab + solvblok->fcblknm;
        SolverCblk *lcblk = solvmtx->cblktab + solvblok->lcblknm;
        if ( fcblk->cblktype & CBLK_IN_LAST ) {
            if ( lcblk->cblktype & CBLK_IN_LAST ) {
                solvblok->inlast = 2;
            }
            else {
                solvblok->inlast = 1;
            }
        }
    }
#else
    (void)solvmtx;
#endif
}

/**
 *******************************************************************************
 *
 * @brief Init a recv solver cblk
 *
 *******************************************************************************
 *
 * @param[in] symbmtx
 *          The pointer to the symbol matrix.
 *
 * @param[inout] solvcblk
 *          The pointer to the cblk.
 *
 * @param[inout] solvblok
 *          The pointer to the first block of the cblk.
 *
 * @param[in] candcblk
 *          Allow us to know the type of the cblk.
 *
 * @param[in] cblklocalnum
 *          Pointer to the cblklocalnum.
 *
 * @param[in] recvidx
 *          Which CBLK_RECV is concerned.
 *
 * @param[in] fcolnum
 *          First column of the cblk.
 *
 * @param[in] lcolnum
 *          Last column of the cblk.
 *
 * @param[in] brownum
 *          Index in th browtab.
 *
 * @param[in] stride
 *          Stride of the cblk.
 *
 * @param[in] cblknum
 *          Global cblk index.
 *
 * @param[in] ownerid
 *          Owner of the cblk.
 *
 *******************************************************************************/
void
solvMatGen_init_cblk_recv( const symbol_matrix_t *symbmtx,
                           const SimuCtrl        *simuctrl,
                                 SolverCblk      *solvcblk,
                                 SolverBlok      *solvblok,
                                 Cand            *candcblk,
                                 pastix_int_t    *cblklocalnum,
                                 pastix_int_t     recvidx,
                                 pastix_int_t     brownum,
                                 pastix_int_t     cblknum,
                                 int              ownerid )
{
    assert( solvblok != NULL );
    assert( brownum >= 0 );

    symbol_cblk_t *symbcblk = symbmtx->cblktab + cblknum;
    symbol_blok_t *symbblok = symbmtx->bloktab + symbcblk->bloknum;

    pastix_int_t lcblknm;
    pastix_int_t stride;
    pastix_int_t layout2D  = candcblk->cblktype & CBLK_LAYOUT_2D;

    assert( symbblok->lcblknm == cblknum );

    lcblknm = cblklocalnum[symbblok->lcblknm] - recvidx;
    stride  = solvMatGen_remote_bloks( symbmtx, symbcblk, simuctrl,
                                       solvcblk, solvblok,
                                       lcblknm, layout2D, ownerid );
    solvblok->fcblknm = cblklocalnum[symbblok->fcblknm];

    solvMatGen_init_cblk( solvcblk, solvblok, candcblk, symbcblk,
                          solvcblk->fcolnum, solvcblk->lcolnum, brownum, stride,
                          cblknum, ownerid );

    solvcblk->brown2d   = brownum;
    solvcblk->cblktype |= CBLK_RECV;

    /* No low-rank compression in distributed for the moment */
    if( solvcblk->cblktype & CBLK_COMPRESSED ) {
        solvcblk->cblktype &= (~CBLK_COMPRESSED);
    }

    /* No Schur complement in distributed for the moment */
    if( solvcblk->cblktype & CBLK_IN_SCHUR ) {
        solvcblk->cblktype &= (~CBLK_IN_SCHUR);
    }
}

/**
 *@}
 */
