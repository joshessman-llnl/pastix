/**
 *
 * @file bcsc_zinit.c
 *
 * @copyright 2004-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.1.0
 * @author Mathieu Faverge
 * @author Pierre Ramet
 * @author Xavier Lacoste
 * @date 2020-01-26
 *
 * @precisions normal z -> c d s
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include "spm.h"
#include "solver.h"
#include "bcsc.h"
#include "bcsc_z.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the values in the block csc stored in the given spm.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc values field is updated.
 *
 *******************************************************************************/
static inline void
bcsc_zinit_A( const spmatrix_t     *spm,
              const pastix_order_t *ord,
              const SolverMatrix   *solvmtx,
              const pastix_int_t   *col2cblk,
                    pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
    pastix_int_t itercblk, iterbcsc, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = spm->dof;

    baseval = spm->colptr[0];

    /**
     * Initialize the values of the matrix A in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t *coltab;
        pastix_int_t  fcolnum, frow, lrow;
        pastix_int_t  itercol2 = ord->permtab[itercol] * dof;
        itercblk = col2cblk[ itercol2 ];

        /* The block column is not stored locally, we skip it */
        if ( itercblk == -1 ) {
            continue;
        }

        fcolnum  = solvmtx->cblktab[itercblk].fcolnum;
        iterbcsc = solvmtx->cblktab[itercblk].bcscnum;
        coltab   = bcsc->cscftab[iterbcsc].coltab;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;
            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol - fcolnum;
                pastix_int_t rowidx = iterrow2;
                pastix_int_t pos = coltab[ colidx ];

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    bcsc->rowtab[ pos ] = rowidx;
                    Lvalues[ pos ] = values[ ival ];
                }

                coltab[ colidx ] += dof;
                assert( coltab[ colidx ] <= coltab[ colidx+1 ] );
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the values in the block csc (upper part) for a symmetric
 * matrix since only one side has been initialized by bcsc_zinit_A()
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc values field is updated.
 *
 *******************************************************************************/
static inline void
bcsc_zinit_Lt( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
                     pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
    pastix_int_t itercblk, iterbcsc, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = spm->dof;

    baseval = spm->colptr[0];

    /**
     * Initialize the values of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ( (itercblk == -1) || (iterrow == itercol) ) {
                continue;
            }

            fcolnum  = solvmtx->cblktab[itercblk].fcolnum;
            iterbcsc = solvmtx->cblktab[itercblk].bcscnum;
            coltab   = bcsc->cscftab[iterbcsc].coltab;

            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    bcsc->rowtab[ pos ] = colidx;
                    Lvalues[ pos ] = values[ ival ];

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the values in the block csc (upper part) for an hermitian
 * matrix since only one side has been initialized by bcsc_zinit_A()
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc values field is updated.
 *
 *******************************************************************************/
#if defined(PRECISION_z) || defined(PRECISION_c)
static inline void
bcsc_zinit_Lh( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
                     pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Lvalues = (pastix_complex64_t*)(bcsc->Lvalues);
    pastix_int_t itercblk, iterbcsc, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = spm->dof;

    baseval = spm->colptr[0];

    /**
     * Initialize the values of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ( (itercblk == -1) || (iterrow == itercol) ) {
                continue;
            }

            fcolnum  = solvmtx->cblktab[itercblk].fcolnum;
            iterbcsc = solvmtx->cblktab[itercblk].bcscnum;
            coltab   = bcsc->cscftab[iterbcsc].coltab;

            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++, pos++)
                {
                    pos = coltab[ rowidx ];

                    bcsc->rowtab[ pos ] = colidx;
                    Lvalues[ pos ] = conj( values[ ival ] );

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

/**
 *******************************************************************************
 *
 * @brief Initialize a value array with the transpose of A that will be used to
 * initialize the coeftab arrays.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[out] trowtab
 *          The row tab associated to the transposition of A.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc Uvalues field is updated.
 *
 *******************************************************************************/
void
bcsc_zinit_At( const spmatrix_t     *spm,
               const pastix_order_t *ord,
               const SolverMatrix   *solvmtx,
               const pastix_int_t   *col2cblk,
                     pastix_int_t   *trowtab,
                     pastix_bcsc_t  *bcsc )
{
    pastix_complex64_t *values  = (pastix_complex64_t*)(spm->values);
    pastix_complex64_t *Uvalues = (pastix_complex64_t*)(bcsc->Uvalues);
    pastix_int_t itercblk, iterbcsc, itercol, baseval;
    pastix_int_t i, ival, idofcol, idofrow;
    int dof = spm->dof;

    baseval = spm->colptr[0];

    /**
     * Initialize the values of the matrix A^t in the blocked csc format. This
     * applies the permutation to the values array.
     */
    for (itercol=0; itercol<spm->gN; itercol++)
    {
        pastix_int_t frow, lrow;
        pastix_int_t itercol2 = ord->permtab[itercol] * dof;

        frow = spm->colptr[itercol]   - baseval;
        lrow = spm->colptr[itercol+1] - baseval;

        for (i=frow; i<lrow; i++)
        {
            pastix_int_t *coltab;
            pastix_int_t fcolnum;
            pastix_int_t iterrow  = spm->rowptr[i]-baseval;
            pastix_int_t iterrow2 = ord->permtab[iterrow] * dof;

            itercblk = col2cblk[ iterrow2 ];

            /* The block column is not stored locally, we skip it */
            if ( itercblk == -1 ) {
                continue;
            }

            fcolnum  = solvmtx->cblktab[itercblk].fcolnum;
            iterbcsc = solvmtx->cblktab[itercblk].bcscnum;
            coltab   = bcsc->cscftab[iterbcsc].coltab;

            ival = i * dof * dof;

            for (idofcol = 0; idofcol < dof; idofcol++)
            {
                pastix_int_t colidx = itercol2 + idofcol;
                pastix_int_t rowidx = iterrow2 - fcolnum;
                pastix_int_t pos;

                for (idofrow = 0; idofrow < dof;
                     idofrow++, ival++, rowidx++)
                {
                    pos = coltab[ rowidx ];

                    trowtab[ pos ] = colidx;
                    Uvalues[ pos ] = values[ ival ];

                    coltab[ rowidx ]++;
                }
            }
        }
    }
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Sort the block csc subarray associated to each column block
 *
 *******************************************************************************
 *
 * @param[in] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *
 * @param[in] rowtab
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] valtab
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 *******************************************************************************/
static inline void
bcsc_zsort( const pastix_bcsc_t *bcsc,
            pastix_int_t        *rowtab,
            pastix_complex64_t  *valtab )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t itercblk, itercol, size;
    void *sortptr[2];

    blockcol = bcsc->cscftab;
    for (itercblk=0; itercblk<bcsc->cscfnbr; itercblk++, blockcol++)
    {
        for (itercol=0; itercol<blockcol->colnbr; itercol++)
        {
            int i;
            sortptr[0] = (void*)(rowtab + blockcol->coltab[itercol]);
            sortptr[1] = (void*)(valtab + blockcol->coltab[itercol]);

            size = blockcol->coltab[itercol+1] - blockcol->coltab[itercol];
            for (i=0; i<size; i++) {
                assert( rowtab[ blockcol->coltab[itercol] + i ] != -1);
            }

            z_qsortIntFloatAsc( sortptr, size );
        }
    }
}

/**
 *******************************************************************************
 *
 * @brief Initialize a centralize pastix_complex64_t block csc.
 *
 *******************************************************************************
 *
 * @param[in] spm
 *          The initial sparse matrix in the spm format.
 *
 * @param[in] ord
 *          The ordering that needs to be applied on the spm to generate the
 *          block csc.
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] col2cblk
 *          Array of matching column with cblk indexes.
 *
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************/
void
bcsc_zinit_centralized( const spmatrix_t     *spm,
                        const pastix_order_t *ord,
                        const SolverMatrix   *solvmtx,
                        const pastix_int_t   *col2cblk,
                              int             initAt,
                              pastix_bcsc_t  *bcsc )
{
    pastix_int_t valuesize;

    bcsc->flttype = spm->flttype;
    valuesize = bcsc_init_centralized_coltab( spm, ord, solvmtx, bcsc );

    /**
     * Initialize the blocked structure of the matrix A
     */
    bcsc_zinit_A( spm, ord, solvmtx, col2cblk, bcsc );
    if ( spm->mtxtype == SpmSymmetric ) {
        bcsc_zinit_Lt( spm, ord, solvmtx, col2cblk, bcsc );
    }
#if defined(PRECISION_z) || defined(PRECISION_c)
    else if ( spm->mtxtype == SpmHermitian ) {
        bcsc_zinit_Lh( spm, ord, solvmtx, col2cblk, bcsc );
    }
#endif /* defined(PRECISION_z) || defined(PRECISION_c) */

    /* Restore the correct coltab arrays */
    bcsc_restore_coltab( bcsc );

    /* Sort the csc */
    bcsc_zsort( bcsc, bcsc->rowtab, bcsc->Lvalues );

    if ( spm->mtxtype == SpmGeneral ) {
	/* A^t is not required if only refinement is performed */
        if (initAt) {
            pastix_int_t *trowtab, i;
            MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char );
            MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

            for (i=0; i<valuesize; i++) {
                trowtab[i] = -1;
            }

            bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

            /* Restore the correct coltab arrays */
            bcsc_restore_coltab( bcsc );

	    /* Sort the transposed csc */
	    bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
	    memFree( trowtab );
        }
    }
    else {
        /* In case of SpmHermitian, conj is applied when used to save memory space */
        bcsc->Uvalues = bcsc->Lvalues;
    }
}

/******************************************************************************
 * Macro: SET_CSC_ROW_VAL(itercblk, therow, thecol, val)          *
 ******************************************************************************
 *                        *
 * Fill next CSC_ROW and CSC_VAL.               *
 *                        *
 * Parameters:                      *
 *   itercblk - Column block index.               *
 *   therow   - Row index.                  *
 *   thecol   - Column index.                 *
 *   val      - Value array.                  *
 *                        *
 ******************************************************************************/
#define SET_CSC_ROW_VAL(itercblk, therow, thecol, val)        \
  _set_csc_row_val(solvmtx, thecsc, itercblk, therow, thecol, \
                   dof, iterdofcol, iterdofrow, iter,         \
                   colidx, strdcol, val)
static inline
void _set_csc_row_val(const SolverMatrix *solvmtx,
                      pastix_bcsc_t          *thecsc,
                      pastix_int_t          itercblk,
                      pastix_int_t          therow,
                      pastix_int_t          thecol,
                      pastix_int_t          dof,
                      pastix_int_t          iterdofcol,
                      pastix_int_t          iterdofrow,
                      pastix_int_t          iter,
                      pastix_int_t          colidx,
                      pastix_int_t          strdcol,
                      double       *val) {
  pastix_int_t fcolnum = solvmtx->cblktab[itercblk].fcolnum;
  pastix_int_t validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow;

  colidx = thecsc->cscftab[itercblk].coltab[thecol-fcolnum]++;

  if (strdcol <= colidx)
    {
      errorPrint("%s:%d colidx %ld >= strdcol %ld",
                 __FILE__, __LINE__, (long)colidx, (long)strdcol);

      EXIT(MOD_SOPALIN, UNKNOWN_ERR);
    }

  thecsc->rowtab[colidx] = therow;
  ((double*)(thecsc->Lvalues))[colidx]  = val[validx];
  thecsc->cscftab[itercblk].coltab[thecol-fcolnum]++;

}

#define SET_CSC_ROW_VAL_CONJ(itercblk, therow, thecol, val) do {    \
    pastix_int_t fcolnum = solvmtx->cblktab[itercblk].fcolnum;               \
    pastix_int_t validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow;   \
                                                                    \
    colidx = thecsc->cscftab[itercblk].coltab[thecol-fcolnum];      \
                                                                    \
    if (strdcol <= colidx)                                          \
    {                                                               \
      errorPrint("%s:%d colidx %ld >= strdcol %ld",                 \
                 __FILE__, __LINE__, (long)colidx, (long)strdcol);  \
                                                                    \
      EXIT(MOD_SOPALIN, UNKNOWN_ERR);                               \
    }                                                               \
                                                                    \
    thecsc->rowtab[colidx] = therow;                               \
    ((double*)(thecsc->Lvalues))[colidx] = val[validx];                 \
    thecsc->cscftab[itercblk].coltab[thecol-fcolnum]++;             \
  } while (0)

/******************************************************************************
 * Macro: SET_TRANS_ROW_VAL(itercblk, therow, thecol, val)          *
 ******************************************************************************
 *                        *
 * Fill next entries for transpose CSC.             *
 *                        *
 * Parameters:                      *
 *   itercblk - Column block index.               *
 *   therow   - Row index.                  *
 *   thecol   - Column index.                 *
 *   val      - Value array.                  *
 *                        *
 ******************************************************************************/

#define SET_TRANS_ROW_VAL(itercblk, therow, thecol, val)        \
_set_trans_row_val(solvmtx, itercblk, therow, thecol,           \
                dof, iterdofcol, iterdofrow, iter,              \
                val, trscltb, trowtab, bcsc)
  
static inline
void _set_trans_row_val(const SolverMatrix* solvmtx, const pastix_int_t itercblk, const pastix_int_t therow, const pastix_int_t thecol, 
                        const pastix_int_t dof, const pastix_int_t iterdofcol, const pastix_int_t iterdofrow, const pastix_int_t iter,
                        const double* val, pastix_int_t** trscltb, pastix_int_t* trowtab, pastix_bcsc_t  *bcsc)
{
    pastix_int_t fcolnum = solvmtx->cblktab[itercblk].fcolnum;
    pastix_int_t trsidx  = trscltb[itercblk][therow-fcolnum];
    pastix_int_t validx  = (iter-1)*dof*dof + iterdofcol*dof + iterdofrow;

    ((double*)(bcsc->Uvalues))[trsidx] = val[validx];
    trowtab[trsidx] = thecol;

    assert(therow >= fcolnum);
    trscltb[itercblk][therow -fcolnum]++;
}
/******************************************************************************
 * Function: CscdOrdistrib                                                    *
 ******************************************************************************
 *                                                                            *
 * Fill in *thecsc* CSC matrix in column block representation.                *
 *                                                                            *
 * - Construct cachetab (sizeof(pastix_int_t)*globalNbCol) which will contain        *
 *   the column block wich will own each column (internal numerotation),      *
 *   or -1 if not local                   *
 *                        *
 * - Build newcoltab (sizeof(pastix_int_t)*globalNbCol) which will contain the         *
 *   coltab corresponding to the local internal CSCd.             *
 *   This CSCd correspond to the given CSCd adding upper part in        *
 *   Symmetric matrix.                    *
 *   Also count number of triples (i,j,v) to send to each other processors.   *
 *                        *
 * - Send the information about how many triples will be sent           *
 *                        *
 * - Fill-in the arrays containing triples to send and send them.         *
 *                        *
 * - Receive those arrays and correct the newcoltab arrays with information   *
 *   from others processors.                  *
 *                        *
 * - Build CSC_COLNBR from symbolic matrix informations and CSC_COL from      *
 *   newcoltab.                     *
 *                        *
 * - Construct transpose matrix, in symmetric mode, transcsc == CSC_VALTAB;   *
 *   in unsymmetric mode, allocate trowtab (number of total local elements),  *
 *   and build trscltb which contains number of elements,           *
 *   in each column of each column bloc.              *
 *                        *
 * - fill-in internal CSC row and values from local given CSCd,         *
 *   also fill-in trowtab and transcsc in unsymmetric mode.           *
 *   CSC_COL and trscltb are incremented for each element added.        *
 *                        *
 * - fill-in  internal CSC row and values from iniformation received,         *
 *   also fill in transposed CSCd in unsymmetric mode.            *
 *   CSC_COL and trscltb are incremented for each element added.        *
 *                        *
 * - restore CSC_COL.                     *
 *                        *
 * - sort internal CSCd.                  *
 *                        *
 * - sort intranal transposed CSCd.                 *
 *                        *
 * Parameters:                      *
 *                                                                            *
 *   thecsc     - Matrix in block column CSC format to fill in.         *
 *   Type       - 3 charactères for matrix Type : only Type[1] is used to     *
 *                check if matrix is Symetric(S) or not(U).           *
 *   transcsc   - Transpose of the CSC in non symetric mode.          *
 *   ord        - ordering                  *
 *   Ncol       - Number of columns.                *
 *   colptr     - Index in *rowind* and *val* of the start of each column.    *
 *   rowind     - Index of the elements.              *
 *   val        - values of the elements.               *
 *   l2g        - global numbers of local nodes.            *
 *   gNcol      - global number of columns.               *
 *   g2l        - local numbers of global nodes, if not local contains -owner *
 *   forcetrans - If matrix symetric, transcsc will be the copy of the        *
 *                CSC_VALTAB.                   *
 *   solvmtx    - Solver matrix                 *
 *   procnum    - MPI process number.                 *
 *   dof        - Number of degree of freedom.              *
 *   comm       - MPI communicator.                 *
 *                                                                            *
 ******************************************************************************/
// void CscdOrdistrib(CscMatrix          *thecsc,
//                    char               *Type,
//                    double             **transcsc,
//                    const Order        *ord,
//                    pastix_int_t                 Ncol,
//                    pastix_int_t                *colptr,
//                    pastix_int_t                *rowind,
//                    double              *val,
//                    pastix_int_t                *l2g,
//                    pastix_int_t                 gNcol,
//                    pastix_int_t                *g2l,
//                    pastix_int_t                 forcetrans,
//                    const SolverMatrix *solvmtx,
//                    pastix_int_t                 procnum,
//                    pastix_int_t                 dof,
//                    MPI_Comm            comm)
void
bcsc_zinit_dist( const spmatrix_t     *spm,
                        const pastix_order_t *ord,
                        const SolverMatrix   *solvmtx,
                        const pastix_int_t   *col2cblk,
                              int             initAt,
                              pastix_bcsc_t  *bcsc )
{
//     pastix_int_t valuesize;

//     bcsc->flttype = spm->flttype;
//     valuesize = bcsc_init_centralized_coltab( spm, ord, solvmtx, bcsc );

//     /**
//      * Initialize the blocked structure of the matrix A
//      */
//     bcsc_zinit_A( spm, ord, solvmtx, col2cblk, bcsc );
//     if ( spm->mtxtype == SpmSymmetric ) {
//         bcsc_zinit_Lt( spm, ord, solvmtx, col2cblk, bcsc );
//     }
// #if defined(PRECISION_z) || defined(PRECISION_c)
//     else if ( spm->mtxtype == SpmHermitian ) {
//         bcsc_zinit_Lh( spm, ord, solvmtx, col2cblk, bcsc );
//     }
// #endif /* defined(PRECISION_z) || defined(PRECISION_c) */

//     /* Restore the correct coltab arrays */
//     bcsc_restore_coltab( bcsc );

//     /* Sort the csc */
//     bcsc_zsort( bcsc, bcsc->rowtab, bcsc->Lvalues );

//     if ( spm->mtxtype == SpmGeneral ) {
// 	/* A^t is not required if only refinement is performed */
//         if (initAt) {
//             pastix_int_t *trowtab, i;
//             MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char );
//             MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

//             for (i=0; i<valuesize; i++) {
//                 trowtab[i] = -1;
//             }

//             bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

//             /* Restore the correct coltab arrays */
//             bcsc_restore_coltab( bcsc );

// 	    /* Sort the transposed csc */
// 	    bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
// 	    memFree( trowtab );
//         }
//     }
//     else {
//         /* In case of SpmHermitian, conj is applied when used to save memory space */
//         bcsc->Uvalues = bcsc->Lvalues;
//     }
    // Old parameters
    pastix_bcsc_t* thecsc = bcsc;
    // Type is handled differently
    // FIXME this is probably wrong
    // transcsc is handled as a member (like Type)
    // ord is already there
    pastix_int_t Ncol = spm->n;
    pastix_int_t* colptr = spm->colptr;
    pastix_int_t* rowind = spm->rowptr;
    double* val = spm->values;
    pastix_int_t* l2g = spm->loc2glob;
    pastix_int_t gNcol = spm->gN;
    pastix_int_t* g2l = spm->glob2loc;
    // solvmtx already there
    pastix_int_t forcetrans = initAt; // This is probably wrong
    pastix_int_t procnum = spm->clustnum;
    pastix_int_t dof = spm->dof;
    SPM_Comm comm = spm->comm;
    
    pastix_int_t          index;
    pastix_int_t          itercol;
    pastix_int_t          newcol;
    // pastix_int_t          loccol;
    pastix_int_t          iter;
    pastix_int_t          rowp1;
    pastix_int_t          colidx;
    pastix_int_t          itercblk;
    pastix_int_t          itercblk2;
    pastix_int_t          strdcol     = 0;
    pastix_int_t        **trscltb     = NULL;
    pastix_int_t         *trowtab     = NULL;
    pastix_int_t         *cachetab    = NULL;
    int          commSize;
    pastix_int_t          proc;
    pastix_int_t         *tosend      = NULL;
    // MPI_Request *tosend_req  = NULL;
    MPI_Request *tosend_creq = NULL;
    MPI_Request *tosend_rreq = NULL;
    MPI_Request *tosend_vreq = NULL;
    pastix_int_t         *torecv      = NULL;
    // MPI_Request *torecv_req  = NULL;
    pastix_int_t         *tosend_cnt  = NULL;
    pastix_int_t        **tosend_col  = NULL;
    pastix_int_t        **tosend_row  = NULL;
    double      **tosend_val  = NULL;
    pastix_int_t        **torecv_col  = NULL;
    pastix_int_t        **torecv_row  = NULL;
    double      **torecv_val  = NULL;
    pastix_int_t         *newcoltab   = NULL;
    // pastix_int_t          owner;
    pastix_int_t          therow;
    #ifndef FORCE_NOMPI
    MPI_Status   status;
    #endif
    pastix_int_t          iterdofrow;
    pastix_int_t          iterdofcol;
    pastix_int_t          nodeidx;
    pastix_int_t          colsize;
    (void)comm;

  
    #ifdef CSC_LOG
    fprintf(stdout, "-> CscdOrdistrib \n");
    #endif
    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gN;
    bcsc->n       = spm->n;
    #if (DBG_SOPALIN_TIME==1)
    Clock clk;
    clockInit(&clk);
    clockStart(&clk);
    #endif

    // FIXME this is almost definitely wrong
    strdcol = bcsc_init_dist_coltab(spm, ord, solvmtx, bcsc);


    /* Initialize some MPI structures */
    MPI_Comm_size(comm, &commSize);
    
    /* tosend will count the number of coefficient to send */
    MALLOC_INTERN(tosend, commSize, pastix_int_t);
    for (proc = 0; proc < commSize; proc++)
        tosend[proc] = 0;
    // MALLOC_INTERN(tosend_req,  commSize, MPI_Request);
    MALLOC_INTERN(tosend_creq, commSize, MPI_Request);
    MALLOC_INTERN(tosend_rreq, commSize, MPI_Request);
    MALLOC_INTERN(tosend_vreq, commSize, MPI_Request);

    /* cachetab: contain the column block or -1 if not local */
    MALLOC_INTERN(cachetab, (gNcol+1)*dof, pastix_int_t);
    for (itercol=0; itercol< (gNcol+1)*dof; itercol++)
        cachetab[itercol] = -1;
    for (itercblk=0; itercblk<solvmtx->cblknbr; itercblk++)
    {
        for (itercol=solvmtx->cblktab[itercblk].fcolnum;
            itercol<solvmtx->cblktab[itercblk].lcolnum+1;
            itercol++)
        {
        cachetab[itercol] = itercblk;
        }
    }

    //   /* newcoltab will contain the number of element in each column
    //    of the symetrized matrix in the new ordering */
    //   MALLOC_INTERN(newcoltab, gNcol+1, pastix_int_t);

    //   for(index=0; index<(gNcol+1); index++)
    //     newcoltab[index] = 0;

    //   for (itercol=0; itercol<Ncol; itercol++)
    //   {
    //     newcol = ord->permtab[l2g[itercol]-1];

    //     newcoltab[newcol] += colptr[itercol+1] - colptr[itercol];

    //     for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
    //     {

    //       loccol = g2l[rowind[iter-1]-1];
    //       if (loccol > 0)
    //       {
    //         if (spm->mtxtype == SpmSymmetric || spm->mtxtype == SpmHermitian)
    //         {
    //           if (loccol -1 != itercol)
    //           {
    //             newcol = ord->permtab[rowind[iter-1]-1];
    //             newcoltab[newcol]++;
    //           }
    //         }
    //       }
    //       else
    //       {
    //         tosend[-loccol]++;
    //       }
    //     }
    //   }

    //   /* Will recv information about what will be received from other processors */
    //   MALLOC_INTERN(torecv, commSize, pastix_int_t);
    //   for (proc = 0; proc < commSize; proc++)
    //     torecv[proc] = 0;

    //   MALLOC_INTERN(torecv_req, commSize, MPI_Request);

    //   for (proc = 0; proc < commSize; proc++)
    //   {
    //     if (proc != procnum)
    //     {
    //       MPI_Isend(&tosend[proc], 1, PASTIX_MPI_INT, proc, 0,
    //                          comm, &tosend_req[proc]);
    //             MPI_Irecv(&torecv[proc], 1, PASTIX_MPI_INT, proc, 0,
    //                          comm, &torecv_req[proc]);
    //           }
    //   }

    //   /* Will contains values from other processors to add to the local
    //    internal CSCD */
    //   MALLOC_INTERN(tosend_col, commSize, pastix_int_t*);
    //   MALLOC_INTERN(tosend_row, commSize, pastix_int_t*);
    //   MALLOC_INTERN(tosend_val, commSize, double*);
    //   MALLOC_INTERN(tosend_cnt, commSize, pastix_int_t);

    //   for (proc = 0; proc < commSize; proc++)
    //   {
    //     if (proc != procnum && tosend[proc] > 0)
    //     {
    //       MALLOC_INTERN(tosend_col[proc], tosend[proc], pastix_int_t);
    //       MALLOC_INTERN(tosend_row[proc], tosend[proc], pastix_int_t);
    //       MALLOC_INTERN(tosend_val[proc], tosend[proc]*dof*dof, double);
    //     }
    //     tosend_cnt[proc] = 0;
    //   }


    //   /* Filling in sending tabs with values and rows*/
    //   for (itercol=0; itercol<Ncol; itercol++)
    //   {
    //     itercblk = cachetab[(ord->permtab[l2g[itercol]-1])*dof];

    //     if (itercblk != -1)
    //     {
    //       /* ok put the value */
    //       for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
    //       {

    //         rowp1 = rowind[iter-1]-1;
    //         newcol = ord->permtab[l2g[itercol]-1];
    //         therow = ord->permtab[rowp1];

    //         itercblk2 = cachetab[therow*dof];

    //         if (itercblk2 != -1)
    //         {

    //         }
    //         else
    //         {
    //           /* Prepare to send row to the owner */
    //           owner = -g2l[ord->peritab[therow]];

    //           tosend_col[owner][tosend_cnt[owner]] = newcol;
    //           tosend_row[owner][tosend_cnt[owner]] = therow;

    //           memcpy(&(tosend_val[owner][tosend_cnt[owner]*dof*dof]),
    //                  &(val[(iter-1)*dof*dof]),
    //                  sizeof(double)*dof*dof);
    //           tosend_cnt[owner]++;
    //         }
    //       }

    //     }
    //     else
    //     {
    //       /* Impossible*/

    //     }
    //   }

    //   /* Sending values to other processors in IJV format. */
    //   for (proc = 0; proc < commSize; proc++)
    //   {
    //     if (proc != procnum)
    //     {
    //       MPI_Wait(&tosend_req[proc], &status);
    //             if (tosend_cnt[proc] > 0)
    //       {
    //         MPI_Isend(tosend_col[proc], (int)(tosend_cnt[proc]),
    //                            PASTIX_MPI_INT, proc, 1, comm, &tosend_creq[proc]);
    //                 MPI_Isend(tosend_row[proc], (int)(tosend_cnt[proc]),
    //                            PASTIX_MPI_INT, proc, 2, comm, &tosend_rreq[proc]);
    //                 MPI_Isend(tosend_val[proc],
    //                            (int)(tosend_cnt[proc]*dof*dof),
    //                            PASTIX_MPI_DOUBLE, proc, 3, comm, &tosend_vreq[proc]);
    //               }
    //       MPI_Wait(&torecv_req[proc], &status);
    //           }
    //   }

    //   memFree_null(tosend_req);
    //   memFree_null(torecv_req);

    //   /* Receiving information from other processors and updating newcoltab */
    //   MALLOC_INTERN(torecv_col, commSize, pastix_int_t*);
    //   MALLOC_INTERN(torecv_row, commSize, pastix_int_t*);
    //   MALLOC_INTERN(torecv_val, commSize, double*);
    //   for (proc = 0; proc < commSize; proc++)
    //   {
    //     if (proc != procnum && torecv[proc] > 0 )
    //     {
    //       MALLOC_INTERN(torecv_col[proc], torecv[proc], pastix_int_t);
    //       MALLOC_INTERN(torecv_row[proc], torecv[proc], pastix_int_t);
    //       MALLOC_INTERN(torecv_val[proc], torecv[proc]*dof*dof, double);
    //       MPI_Recv(torecv_col[proc], torecv[proc], PASTIX_MPI_INT,
    //                         proc, 1, comm, &status );
    //             MPI_Recv(torecv_row[proc], torecv[proc], PASTIX_MPI_INT,
    //                         proc, 2, comm, &status );
    //             MPI_Recv(torecv_val[proc], torecv[proc]*dof*dof, PASTIX_MPI_DOUBLE,
    //                         proc, 3, comm, &status );
        
    //       if (spm->mtxtype == SpmSymmetric || spm->mtxtype == SpmHermitian)
    //       {
    //         for (iter = 0; iter < torecv[proc]; iter++)
    //         {
    //           newcol= torecv_row[proc][iter];
    //           newcoltab[newcol] ++;
    //         }
    //       }
    //     }
    //   }


    // #if (DBG_SOPALIN_TIME==1)
    //   clockStop(&(clk));
    //   fprintf(stdout, "CscdOrdistrib step 1 : %.3g s\n",
    //           (double)clockVal(&clk));
    //   clockInit(&clk);
    //   clockStart(&clk);
    // #endif
    //   /* Finishing newcoltab construction :
    //    *
    //    * Now, newcoltab will contain starting index of each
    //    * column of rows and values in new ordering
    //    */
    //   newcol = 0;
    //   for (index=0; index<(gNcol+1); index++)
    //   {
    //     colidx = newcoltab[index];
    //     newcoltab[index] = newcol;
    //     newcol += colidx;
    //   }

    // SET_CSC_COL(newcoltab);
    /* Local coltab */
        // thecsc->cscfnbr = solvmtx->cblknbr;
        // MALLOC_INTERN(thecsc->cscftab, thecsc->cscfnbr, bcsc_cblk_t);

        // for (index=0; index<solvmtx->cblknbr; index++)
        // {
        //   pastix_int_t fcolnum = solvmtx->cblktab[index].fcolnum;
        //   pastix_int_t lcolnum = solvmtx->cblktab[index].lcolnum;
        //   thecsc->cscftab[index].colnbr = (lcolnum - fcolnum+1);

        //   MALLOC_INTERN(thecsc->cscftab[index].coltab,
        //                 thecsc->cscftab[index].colnbr + 1, pastix_int_t);

        //   if ((fcolnum)%dof != 0)
        //     errorPrint("dof doesn't divide fcolnum");

        //   colsize = 0;
        //   for (iter=0; iter<(thecsc->cscftab[index].colnbr + 1); iter++)
        //   {
        //     /* fcolnum %dof = 0 */
        //     nodeidx = (fcolnum+(iter-iter%dof))/dof;
        //     if (g2l != NULL &&
        //         iter != thecsc->cscftab[index].colnbr &&
        //         !(g2l[ord->peritab[nodeidx]] > 0))
        //     {
        //       errorPrint("Columns in internal CSCD must be in given CSCD");
        //     }

        //     thecsc->cscftab[index].coltab[iter] = colsize+ strdcol;
        //     if (iter < thecsc->cscftab[index].colnbr)
        //     {
        //       colsize = (newcoltab[nodeidx+1] -
        //                  newcoltab[nodeidx])*dof;
        //     }
        //     strdcol =  thecsc->cscftab[index].coltab[iter];
        //   }

        //   /* strdcol <- colptr[n] for the index_th column block */
        //   /*   ie : the number of element in the column block */

        //   strdcol = thecsc->cscftab[index].coltab[thecsc->cscftab[index].colnbr];
        // }
        // memFree_null(newcoltab);

    #if (DBG_SOPALIN_TIME==1)
    clockStop(&(clk));
    fprintf(stdout, "CscdOrdistrib step 2 : %.3g s\n",
            (double)clockVal(&clk));
    clockInit(&clk);
    clockStart(&clk);
    #endif

    // CSC_ALLOC macro
    // MALLOC_INTERN(thecsc->rowtab, strdcol, pastix_int_t);              
    // MALLOC_INTERN(thecsc->Lvalues, strdcol, double);            // FIXME: Transpose values?
                                                                    
        // if (bcsc->Uvalues != NULL)                                         
        // {     
        if (spm->mtxtype == SpmSymmetric || spm->mtxtype == SpmHermitian)
        {                                                           
            if (forcetrans == 1)                                
            {                                                         
                thecsc->Uvalues = thecsc->Lvalues;                  
            }                                                         
        }                                                           
        else                                                        
        {                                                           
            MALLOC_INTERN( bcsc->Uvalues, strdcol * pastix_size_of( bcsc->flttype ), char );          
            MALLOC_INTERN(trowtab, strdcol, pastix_int_t);                     
            MALLOC_INTERN(trscltb, solvmtx->cblknbr, pastix_int_t *);          
                                                                        
            for (index=0; index<solvmtx->cblknbr; index++)            
            {                                                         
                MALLOC_INTERN(trscltb[index],                           
                            thecsc->cscftab[index].colnbr + 1, pastix_int_t);         
                for (iter=0; iter<(thecsc->cscftab[index].colnbr + 1); iter++) 
                {                                                       
                    trscltb[index][iter] = thecsc->cscftab[index].coltab[iter];    
                }                                                       
            }                                                         
        }                                                           
        // }                                                             

    /* Filling in thecsc with values and rows*/
    for (itercol=0; itercol<Ncol; itercol++)
    {
        itercblk = cachetab[(ord->permtab[l2g[itercol]-1])*dof];

        if (itercblk != -1)
        {
        /* ok put the value */
        for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
        {

            for (iterdofcol = 0; iterdofcol < dof; iterdofcol++)
            {
            for (iterdofrow = 0; iterdofrow < dof; iterdofrow++)
            {
                rowp1 = rowind[iter-1]-1;
                therow = ord->permtab[rowp1]*dof + iterdofrow;
                newcol = (ord->permtab[l2g[itercol]-1])*dof+iterdofcol;
                SET_CSC_ROW_VAL(itercblk, therow, newcol, val);

                itercblk2 = cachetab[therow];

                if (itercblk2 != -1)
                {
                switch (spm->mtxtype)
                {
                case SpmSymmetric:
                {
                    if (rowp1 != l2g[itercol]-1)
                    {
                    /* same thing but newcol <-> therow */
                    SET_CSC_ROW_VAL(itercblk2, newcol, therow,
                                    val);

                    }
                    break;
                }
                #if defined(PRECISION_z) || defined(PRECISION_c)
                case SpmHermitian:
                {
                    if (rowp1 != l2g[itercol]-1)
                    {
                    /* same thing but newcol <-> therow */
                    SET_CSC_ROW_VAL_CONJ(itercblk2, newcol, therow,
                                    val);

                    }
                    break;
                }
                #endif /* #if defined(PRECISION_z) || defined(PRECISION_c) */
                case SpmGeneral:
                {
                    // if (bcsc->Uvalues != NULL)
                    // {
                    SET_TRANS_ROW_VAL(itercblk2, therow, newcol,
                                        val);
                    // }
                    break;
                }
                default:
                    errorPrint("Unknown Matrix Type");
                    EXIT(MOD_SOPALIN, UNKNOWN_ERR);
                }

                }
            }
            }
        }

        }
        else
        {
        /* Impossible*/
        errorPrint("Error in CscdOrdistrib");
        EXIT(MOD_SOPALIN, UNKNOWN_ERR)
            }
    }

    memFree_null(tosend);

    for (proc = 0; proc < commSize; proc++)
    {
        if (proc != procnum)
        {
        for (iter = 0; iter < torecv[proc]; iter++)
        {
            for (iterdofcol = 0; iterdofcol < dof; iterdofcol++)
            {
            for (iterdofrow = 0; iterdofrow < dof; iterdofrow++)
            {
                switch (spm->mtxtype)
                {
                case SpmSymmetric:
                {
                newcol  = torecv_col[proc][iter]*dof+iterdofcol;
                therow  = torecv_row[proc][iter]*dof+iterdofrow;
                itercblk2 = cachetab[therow];
                /* iter is 0 based here, not in SET_CSC_ROW_VAL */
                iter++;
                SET_CSC_ROW_VAL(itercblk2, newcol, therow,
                                torecv_val[proc]);
                iter--;
                break;
                }
                #if defined(PRECISION_z) || defined(PRECISION_c)
                case SpmHermitian:
                {
                newcol  = torecv_col[proc][iter]*dof+iterdofcol;
                therow  = torecv_row[proc][iter]*dof+iterdofrow;
                itercblk2 = cachetab[therow];
                /* iter is 0 based here, not in SET_CSC_ROW_VAL */
                iter++;
                SET_CSC_ROW_VAL_CONJ(itercblk2, newcol, therow,
                                    torecv_val[proc]);
                iter--;
                break;
                }
                #endif /* #if defined(PRECISION_z) || defined(PRECISION_c) */
                case SpmGeneral:
                {
                // if (bcsc->Uvalues != NULL)
                // {
                    newcol = torecv_col[proc][iter]*dof+iterdofcol;
                    therow = torecv_row[proc][iter]*dof+iterdofrow;
                    itercblk2 = cachetab[therow];
                    /* iter is 0 based here,
                    not in SET_TRANS_ROW_VAL */
                    iter++;
                    SET_TRANS_ROW_VAL(itercblk2, therow, newcol,
                                    torecv_val[proc]);
                    iter--;
                // }
                break;
                }
                default:
                errorPrint("Unknown Matrix Type");
                EXIT(MOD_SOPALIN, UNKNOWN_ERR);
                }
            }
            }
        }
        if (torecv[proc] > 0)
        {
            memFree_null(torecv_col[proc]);
            memFree_null(torecv_row[proc]);
            memFree_null(torecv_val[proc]);
        }
        }
    }
    memFree_null(torecv_col);
    memFree_null(torecv_row);
    memFree_null(torecv_val);
    memFree_null(cachetab);
    memFree_null(torecv);
    for (proc = 0; proc < commSize; proc++)
    {
        if (proc != procnum)
        {
        if (tosend_cnt[proc] > 0)
        {
            MPI_Wait(&tosend_creq[proc], &status);
                    memFree_null(tosend_col[proc]);
            MPI_Wait(&tosend_rreq[proc], &status);
                    memFree_null(tosend_row[proc]);
            MPI_Wait(&tosend_vreq[proc], &status);
                    memFree_null(tosend_val[proc]);
        }
        }
    }
    memFree_null(tosend_creq);
    memFree_null(tosend_rreq);
    memFree_null(tosend_vreq);
    memFree_null(tosend_cnt);
    memFree_null(tosend_col);
    memFree_null(tosend_row);
    memFree_null(tosend_val);

    if (trscltb != NULL)
    {
        for (index=0; index<solvmtx->cblknbr; index++)
        {
        memFree_null(trscltb[index]);
        }
        memFree_null(trscltb);
    }

    /* 2nd membre */
    /* restore good coltab */
    colidx = 0;
    for (index=0; index<solvmtx->cblknbr; index++)
    {
        for(iter=0;iter<(thecsc->cscftab[index].colnbr+1); iter++)
        {
        newcol = thecsc->cscftab[index].coltab[iter];
        thecsc->cscftab[index].coltab[iter] = colidx;
        colidx = newcol;
        }
    }

    #if (DBG_SOPALIN_TIME==1)
    clockStop(&(clk));
    fprintf(stdout, "CscdOrdistrib step 3 : %.3g s\n",
            (double)clockVal(&clk));
    clockInit(&clk);
    clockStart(&clk);
    #endif

    //   CSC_SORT;
    /* Sort */                                                         
        // for (index=0; index<solvmtx->cblknbr; index++) {                   
        //   for (iter=0; iter<thecsc->cscftab[index].colnbr; iter++) {            
        //     pastix_int_t   *t = &(thecsc->rowtab[thecsc->cscftab[index].coltab[iter]]);              
        //     double *v = &(((double*)thecsc->Lvalues)[thecsc->cscftab[index].coltab[iter]]);              
        //     pastix_int_t    n = thecsc->cscftab[index].coltab[iter+1]-                
        //       thecsc->cscftab[index].coltab[iter];                                  
        //     pastix_int_t    ndof2 = 1; /* internal CSC is with one DoF */    
        //     void * sortptr[3];                                             
        //     sortptr[0] = t;                                                
        //     sortptr[1] = v;                                                
        //     sortptr[2] = &ndof2;                                           
        //     z_qsortIntFloatAsc(sortptr, n);                                  
        //   }                                                                
        // }         
           
        // if (spm->mtxtype == SpmGeneral && forcetrans) {
        //     bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );
        // }
        // if (bcsc->Uvalues != NULL) {                     
        // if (!forcetrans) {                                               
        //     // for (index=0; index<solvmtx->cblknbr; index++) {               
        //     //   for (iter=0; iter<thecsc->cscftab[index].colnbr; iter++) {        
        //     //     pastix_int_t   *t = &(trowtab[thecsc->cscftab[index].coltab[iter]]);  
        //     //     double *v = &(((double*)bcsc->Uvalues)[thecsc->cscftab[index].coltab[iter]]);     
        //     //     pastix_int_t n;                                              
        //     //     pastix_int_t    ndof2 = 1; /* internal CSC is with 1 DoF */  
        //     //     void * sortptr[3];                                         
                                                                        
        //     //     n = thecsc->cscftab[index].coltab[iter+1] -                         
        //     //       thecsc->cscftab[index].coltab[iter];                              
                                                                        
        //     //     sortptr[0] = t;                                            
        //     //     sortptr[1] = v;                                            
        //     //     sortptr[2] = &ndof2;                                       
        //     //     z_qsortIntFloatAsc(sortptr, n);                              
        //     //   }                                                            
        //     // }                                                              
        //     // memFree_null(trowtab);                
        //     // pastix_int_t *trowtab, i;
        //     MALLOC_INTERN( bcsc->Uvalues, strdcol * pastix_size_of( bcsc->flttype ), char );
        //     // MALLOC_INTERN( trowtab, strdcol, pastix_int_t);

        //     // for (i=0; i<strdcol; i++) {
        //     //     trowtab[i] = -1;
        //     // }

        //     // bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

        //     /* Restore the correct coltab arrays */
        //     bcsc_restore_coltab( bcsc );

        //     /* Sort the transposed csc */
        //     bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
        //     memFree( trowtab );                         
        // }   
        // }     
            /**
     * Initialize the blocked structure of the matrix A
     */
        bcsc_zinit_A( spm, ord, solvmtx, col2cblk, bcsc );
        if ( spm->mtxtype == SpmSymmetric ) {
            bcsc_zinit_Lt( spm, ord, solvmtx, col2cblk, bcsc );
        }
        #if defined(PRECISION_z) || defined(PRECISION_c)
            else if ( spm->mtxtype == SpmHermitian ) {
                bcsc_zinit_Lh( spm, ord, solvmtx, col2cblk, bcsc );
            }
        #endif /* defined(PRECISION_z) || defined(PRECISION_c) */

        /* Restore the correct coltab arrays */
        bcsc_restore_coltab( bcsc );

        bcsc_zsort( bcsc, bcsc->rowtab, bcsc->Lvalues );

        if ( spm->mtxtype == SpmGeneral ) {
        /* A^t is not required if only refinement is performed */
            if (initAt) {
                // pastix_int_t *trowtab, i;
                // MALLOC_INTERN( bcsc->Uvalues, valuesize * pastix_size_of( bcsc->flttype ), char );
                // MALLOC_INTERN( trowtab, valuesize, pastix_int_t);

                // for (i=0; i<valuesize; i++) {
                //     trowtab[i] = -1;
                // }

                bcsc_zinit_At( spm, ord, solvmtx, col2cblk, trowtab, bcsc );

                /* Restore the correct coltab arrays */
                bcsc_restore_coltab( bcsc );

            /* Sort the transposed csc */
            bcsc_zsort( bcsc, trowtab, bcsc->Uvalues );
            memFree( trowtab );
            }
        }
        else {
            /* In case of SpmHermitian, conj is applied when used to save memory space */
            bcsc->Uvalues = bcsc->Lvalues;
        }                                                             

    #if (DBG_SOPALIN_TIME==1)
    clockStop(&(clk));
    fprintf(stdout, "CscdOrdistrib step 4 : %.3g s\n",
            (double)clockVal(&clk));
    #endif
    #ifdef CSC_LOG
    fprintf(stdout, "<- CscdOrdistrib \n");
    #endif
}
