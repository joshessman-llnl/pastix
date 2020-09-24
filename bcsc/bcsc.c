/**
 *
 * @file bcsc.c
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
 **/
#include "common.h"
#include "pastix/order.h"
#include "spm.h"
#include "solver.h"
#include "bcsc.h"

#include "bcsc_z.h"
#include "bcsc_c.h"
#include "bcsc_d.h"
#include "bcsc_s.h"

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Initialize the coltab of a block csc matrix.
 *
 *******************************************************************************
 *
 * @param[in] solvmtx
 *          The solver matrix structure that describe the data distribution.
 *
 * @param[in] newcoltab
 *          Array of size spm->gN+1. This array is global coltab with -1 for non
 *          local indexes.
 *
 * @param[in] dof
 *          The degree of freedom of each unknown.
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the ininitialized coltab split per block.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
static inline pastix_int_t
bcsc_init_coltab( const SolverMatrix  *solvmtx,
                  const pastix_int_t  *newcoltab,
                        pastix_int_t   dof,
                        pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t cblknum, bcscnum, iter, idxcol, nodeidx, colsize;

    SolverCblk *cblk = solvmtx->cblktab;
    bcsc->cscfnbr    = solvmtx->cblknbr - solvmtx->faninnbr - solvmtx->recvnbr;
    MALLOC_INTERN( bcsc->cscftab, bcsc->cscfnbr, bcsc_cblk_t );

    idxcol   = 0;
    cblk     = solvmtx->cblktab;
    blockcol = bcsc->cscftab;
    bcscnum  = 0;
    for (cblknum = 0; cblknum < solvmtx->cblknbr; cblknum++, cblk++)
    {
        pastix_int_t fcolnum = cblk->fcolnum;

        if ( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ) {
            continue;
        }

        blockcol->cblknum = cblknum;
        blockcol->colnbr  = cblk_colnbr( cblk );
        assert( cblk->bcscnum == bcscnum );
        MALLOC_INTERN( blockcol->coltab, blockcol->colnbr + 1, pastix_int_t );

        /* Works only for DoF constant */
        assert( fcolnum % dof == 0 );

        blockcol->coltab[0] = idxcol;
        for (iter=0; iter < blockcol->colnbr; iter++)
        {
            nodeidx = ( fcolnum + (iter-iter%dof) ) / dof;

            colsize = (newcoltab[nodeidx+1] - newcoltab[nodeidx]) * dof;
            blockcol->coltab[iter+1] = blockcol->coltab[iter] + colsize;
        }

        idxcol = blockcol->coltab[blockcol->colnbr];

        blockcol++;
        bcscnum++;
    }
    assert( (blockcol - bcsc->cscftab) == bcsc->cscfnbr );
    assert( bcscnum == bcsc->cscfnbr );

    if ( idxcol > 0 ) {
        MALLOC_INTERN( bcsc->rowtab,  idxcol, pastix_int_t);
        MALLOC_INTERN( bcsc->Lvalues, idxcol * pastix_size_of( bcsc->flttype ), char );
    }
    else {
        bcsc->rowtab  = NULL;
        bcsc->Lvalues = NULL;
    }
    bcsc->Uvalues = NULL;

    return idxcol;
}

/**
 *******************************************************************************
 *
 * @ingroup bcsc_internal
 *
 * @brief Restore the coltab array
 *
 * Function to restore the coltab array when it has been modified to initialize
 * the row and values arrays.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          On entry, the bcsc to restore.
 *          On exit, the coltab array of the bcsc is restored to the correct
 *          indexes.
 *
 *******************************************************************************/
void
bcsc_restore_coltab( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *blockcol;
    pastix_int_t index, iter, idxcol, idxcoltmp;

    idxcol = 0;
    blockcol = bcsc->cscftab;
    for (index=0; index<bcsc->cscfnbr; index++, blockcol++)
    {
        for (iter=0; iter <= blockcol->colnbr; iter++)
        {
            idxcoltmp = blockcol->coltab[iter];
            blockcol->coltab[iter] = idxcol;
            idxcol = idxcoltmp;
        }
    }
    return;
}

/**
 *******************************************************************************
 *
 * @brief Initialize the coltab of a centralized block csc matrix.
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
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
pastix_int_t
bcsc_init_centralized_coltab( const spmatrix_t     *spm,
                              const pastix_order_t *ord,
                              const SolverMatrix   *solvmtx,
                                    pastix_bcsc_t  *bcsc )
{
    pastix_int_t  valuesize, baseval;
    pastix_int_t *globcol  = NULL;
    pastix_int_t *colptr = spm->colptr;
    pastix_int_t *rowptr = spm->rowptr;
    int dof = spm->dof;
    int sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    bcsc->mtxtype = spm->mtxtype;
    baseval = spm->colptr[0];

    /*
     * Allocate and initialize globcol that contains the number of elements in
     * each column of the input matrix
     * Globcol is equivalent to the classic colptr for the internal blocked
     * csc. The blocked csc integrate the perumtation computed within order
     * structure.
     */
    MALLOC_INTERN( globcol, spm->gN+1, pastix_int_t );
    memset( globcol, 0, (spm->gN+1) * sizeof(pastix_int_t) );

    assert( spm->loc2glob == NULL );

    {
        pastix_int_t itercol, newcol;

        for (itercol=0; itercol<spm->gN; itercol++)
        {
            pastix_int_t frow = colptr[itercol]   - baseval;
            pastix_int_t lrow = colptr[itercol+1] - baseval;
            newcol = ord->permtab[itercol];
            globcol[newcol] += lrow - frow;

            assert( (lrow - frow) >= 0 );
            if (sym) {
                pastix_int_t iterrow, newrow;

                for (iterrow=frow; iterrow<lrow; iterrow++)
                {
                    pastix_int_t tmprow = rowptr[iterrow] - baseval;
                    if (tmprow != itercol) {
                        newrow = ord->permtab[tmprow];
                        globcol[newrow]++;
                    }
                }
            }
        }

        /* Compute displacements to update the colptr array */
        {
            pastix_int_t tmp, idx;

            idx = 0;
            for (itercol=0; itercol<=spm->gN; itercol++)
            {
                tmp = globcol[itercol];
                globcol[itercol] = idx;
                idx += tmp;
            }
        }
    }

    valuesize = bcsc_init_coltab( solvmtx, globcol, dof, bcsc );
    memFree_null( globcol );

    return valuesize;
}

/**
 *******************************************************************************
 *
 * @brief Initialize a centralized block csc when no MPI processes are involved.
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
bcsc_init_centralized( const spmatrix_t     *spm,
                       const pastix_order_t *ord,
                       const SolverMatrix   *solvmtx,
                             pastix_int_t    initAt,
                             pastix_bcsc_t  *bcsc )
{
    pastix_int_t  itercol, itercblk;
    pastix_int_t  cblknbr  = solvmtx->cblknbr;
    pastix_int_t  eltnbr   = spm->gNexp;
    pastix_int_t *col2cblk = NULL;

    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gN;
    bcsc->n       = spm->n;

    assert( spm->loc2glob == NULL );

    /*
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, eltnbr, pastix_int_t );
        for (itercol=0; itercol<eltnbr; itercol++)
        {
            col2cblk[itercol] = -1;
        }

        for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
        {
            if( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ){
                continue;
            }
            for (itercol  = cblk->fcolnum;
                 itercol <= cblk->lcolnum;
                 itercol++ )
            {
                col2cblk[itercol] = itercblk;
            }
        }
    }

    /*
     * Fill in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
    switch( spm->flttype ) {
    case SpmFloat:
        bcsc_sinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmDouble:
        bcsc_dinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex32:
        bcsc_cinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex64:
        bcsc_zinit_centralized( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmPattern:
    default:
        fprintf(stderr, "bcsc_init_centralized: Error unknown floating type for input spm\n");
    }

    memFree_null(col2cblk);
}

/**
 *******************************************************************************
 *
 * @brief Initialize the coltab of a centralized block csc matrix.
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
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The number of non zero unknowns in the matrix.
 *
 *******************************************************************************/
pastix_int_t
bcsc_init_dist_coltab( const spmatrix_t     *spm,
                              const pastix_order_t *ord,
                              const SolverMatrix   *solvmtx,
                              const pastix_int_t   *col2cblk,
                                    pastix_bcsc_t  *bcsc )
{
    // pastix_int_t  valuesize, baseval;
    // pastix_int_t *globcol  = NULL;
    // pastix_int_t *colptr = spm->colptr;
    // pastix_int_t *rowptr = spm->rowptr;
    // int dof = spm->dof;
    // int sym = (spm->mtxtype == SpmSymmetric) || (spm->mtxtype == SpmHermitian);

    bcsc->mtxtype = spm->mtxtype;
    // baseval = spm->colptr[0];

    /*
     * Allocate and initialize newcoltab that contains the number of elements in
     * each column of the input matrix
     * Globcol is equivalent to the classic colptr for the internal blocked
     * csc. The blocked csc integrate the perumtation computed within order
     * structure.
     */
    pastix_int_t* colptr = spm->colptr;
    pastix_int_t* rowind = spm->rowptr;
    double* val = spm->values;
    pastix_int_t* l2g = spm->loc2glob;
    pastix_int_t gNcol = spm->gN;
    pastix_int_t Ncol = spm->n;
    pastix_int_t* g2l = spm->glob2loc;
    pastix_int_t procnum = spm->clustnum;
    pastix_int_t dof = spm->dof;
    SPM_Comm comm = spm->comm;
    pastix_int_t  valuesize;
    
    pastix_int_t          index;
    pastix_int_t          itercol;
    pastix_int_t          newcol;
    pastix_int_t          loccol;
    pastix_int_t          iter;
    pastix_int_t          rowp1;
    pastix_int_t          colidx;
    pastix_int_t          itercblk;
    pastix_int_t          itercblk2;
    // pastix_int_t          strdcol     = 0;
    // pastix_int_t        **trscltb     = NULL;
    // pastix_int_t         *trowtab     = NULL;
    int          commSize;
    pastix_int_t          proc;
    pastix_int_t         *tosend      = NULL;
    MPI_Request *tosend_req  = NULL;
    MPI_Request *tosend_creq = NULL;
    MPI_Request *tosend_rreq = NULL;
    MPI_Request *tosend_vreq = NULL;
    pastix_int_t         *torecv      = NULL;
    MPI_Request *torecv_req  = NULL;
    pastix_int_t         *tosend_cnt  = NULL;
    pastix_int_t        **tosend_col  = NULL;
    pastix_int_t        **tosend_row  = NULL;
    double      **tosend_val  = NULL;
    pastix_int_t        **torecv_col  = NULL;
    pastix_int_t        **torecv_row  = NULL;
    double      **torecv_val  = NULL;
    pastix_int_t         *newcoltab   = NULL;
    pastix_int_t          owner;
    pastix_int_t          therow;
    #ifndef FORCE_NOMPI
      MPI_Status   status;
    #endif

    // MPI setup
    MPI_Comm_size(comm, &commSize);
  
    /* tosend will count the number of coefficient to send */
    MALLOC_INTERN(tosend, commSize, pastix_int_t);
    for (proc = 0; proc < commSize; proc++)
      tosend[proc] = 0;
    MALLOC_INTERN(tosend_req,  commSize, MPI_Request);
    MALLOC_INTERN(tosend_creq, commSize, MPI_Request);
    MALLOC_INTERN(tosend_rreq, commSize, MPI_Request);
    MALLOC_INTERN(tosend_vreq, commSize, MPI_Request);

    /* newcoltab will contain the number of element in each column
    of the symetrized matrix in the new ordering */
    MALLOC_INTERN(newcoltab, gNcol+1, pastix_int_t);

    for(index=0; index<(gNcol+1); index++)
      newcoltab[index] = 0;

    for (itercol=0; itercol<Ncol; itercol++)
    {
      newcol = ord->permtab[l2g[itercol]-1];

      newcoltab[newcol] += colptr[itercol+1] - colptr[itercol];

      for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
      {

        loccol = g2l[rowind[iter-1]-1];
        if (loccol > 0)
        {
          if (spm->mtxtype == SpmSymmetric || spm->mtxtype == SpmHermitian)
          {
            if (loccol -1 != itercol)
            {
              newcol = ord->permtab[rowind[iter-1]-1];
              newcoltab[newcol]++;
            }
          }
        }
        else
        {
          tosend[-loccol]++;
        }
      }
    }

    /* Will recv information about what will be received from other processors */
    MALLOC_INTERN(torecv, commSize, pastix_int_t);
    for (proc = 0; proc < commSize; proc++)
      torecv[proc] = 0;

    MALLOC_INTERN(torecv_req, commSize, MPI_Request);

    for (proc = 0; proc < commSize; proc++)
    {
      if (proc != procnum)
      {
        MPI_Isend(&tosend[proc], 1, PASTIX_MPI_INT, proc, 0,
                          comm, &tosend_req[proc]);
              MPI_Irecv(&torecv[proc], 1, PASTIX_MPI_INT, proc, 0,
                          comm, &torecv_req[proc]);
            }
    }

    /* Will contains values from other processors to add to the local
    internal CSCD */
    MALLOC_INTERN(tosend_col, commSize, pastix_int_t*);
    MALLOC_INTERN(tosend_row, commSize, pastix_int_t*);
    MALLOC_INTERN(tosend_val, commSize, double*);
    MALLOC_INTERN(tosend_cnt, commSize, pastix_int_t);

    for (proc = 0; proc < commSize; proc++)
    {
      if (proc != procnum && tosend[proc] > 0)
      {
        MALLOC_INTERN(tosend_col[proc], tosend[proc], pastix_int_t);
        MALLOC_INTERN(tosend_row[proc], tosend[proc], pastix_int_t);
        MALLOC_INTERN(tosend_val[proc], tosend[proc]*dof*dof, double);
      }
      tosend_cnt[proc] = 0;
    }


    /* Filling in sending tabs with values and rows*/
    for (itercol=0; itercol<Ncol; itercol++)
    {
      itercblk = col2cblk[(ord->permtab[l2g[itercol]-1])*dof];

      if (itercblk != -1)
      {
        /* ok put the value */
        for (iter=colptr[itercol]; iter<colptr[itercol+1]; iter++)
        {

          rowp1 = rowind[iter-1]-1;
          newcol = ord->permtab[l2g[itercol]-1];
          therow = ord->permtab[rowp1];

          itercblk2 = col2cblk[therow*dof];

          if (itercblk2 != -1)
          {

          }
          else
          {
            /* Prepare to send row to the owner */
            owner = -g2l[ord->peritab[therow]];

            tosend_col[owner][tosend_cnt[owner]] = newcol;
            tosend_row[owner][tosend_cnt[owner]] = therow;

            memcpy(&(tosend_val[owner][tosend_cnt[owner]*dof*dof]),
                  &(val[(iter-1)*dof*dof]),
                  sizeof(double)*dof*dof);
            tosend_cnt[owner]++;
          }
        }

      }
      else
      {
        /* Impossible*/

      }
    }

    /* Sending values to other processors in IJV format. */
    for (proc = 0; proc < commSize; proc++)
    {
      if (proc != procnum)
      {
        MPI_Wait(&tosend_req[proc], &status);
              if (tosend_cnt[proc] > 0)
        {
          MPI_Isend(tosend_col[proc], (int)(tosend_cnt[proc]),
                            PASTIX_MPI_INT, proc, 1, comm, &tosend_creq[proc]);
                  MPI_Isend(tosend_row[proc], (int)(tosend_cnt[proc]),
                            PASTIX_MPI_INT, proc, 2, comm, &tosend_rreq[proc]);
                  MPI_Isend(tosend_val[proc],
                            (int)(tosend_cnt[proc]*dof*dof),
                            PASTIX_MPI_DOUBLE, proc, 3, comm, &tosend_vreq[proc]);
                }
        MPI_Wait(&torecv_req[proc], &status);
            }
    }

    memFree_null(tosend_req);
    memFree_null(torecv_req);

    /* Receiving information from other processors and updating newcoltab */
    MALLOC_INTERN(torecv_col, commSize, pastix_int_t*);
    MALLOC_INTERN(torecv_row, commSize, pastix_int_t*);
    MALLOC_INTERN(torecv_val, commSize, double*);
    for (proc = 0; proc < commSize; proc++)
    {
      if (proc != procnum && torecv[proc] > 0 )
      {
        MALLOC_INTERN(torecv_col[proc], torecv[proc], pastix_int_t);
        MALLOC_INTERN(torecv_row[proc], torecv[proc], pastix_int_t);
        MALLOC_INTERN(torecv_val[proc], torecv[proc]*dof*dof, double);
        MPI_Recv(torecv_col[proc], torecv[proc], PASTIX_MPI_INT,
                          proc, 1, comm, &status );
              MPI_Recv(torecv_row[proc], torecv[proc], PASTIX_MPI_INT,
                          proc, 2, comm, &status );
              MPI_Recv(torecv_val[proc], torecv[proc]*dof*dof, PASTIX_MPI_DOUBLE,
                          proc, 3, comm, &status );
        
        if (spm->mtxtype == SpmSymmetric || spm->mtxtype == SpmHermitian)
        {
          for (iter = 0; iter < torecv[proc]; iter++)
          {
            newcol= torecv_row[proc][iter];
            newcoltab[newcol] ++;
          }
        }
      }
    }


  #if (DBG_SOPALIN_TIME==1)
    clockStop(&(clk));
    fprintf(stdout, "CscdOrdistrib step 1 : %.3g s\n",
            (double)clockVal(&clk));
    clockInit(&clk);
    clockStart(&clk);
  #endif
    /* Finishing newcoltab construction :
    *
    * Now, newcoltab will contain starting index of each
    * column of rows and values in new ordering
    */
    newcol = 0;
    for (index=0; index<(gNcol+1); index++)
    {
      colidx = newcoltab[index];
      newcoltab[index] = newcol;
      newcol += colidx;
    }

    for (proc = 0; proc < commSize; proc++)
    {
      if (proc != procnum)
      {
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

    valuesize = bcsc_init_coltab( solvmtx, newcoltab, dof, bcsc );
    memFree_null( newcoltab );


    return valuesize;
}



void bcsc_init_dist( const spmatrix_t     *spm,
                       const pastix_order_t *ord,
                       const SolverMatrix   *solvmtx,
                             pastix_int_t    initAt,
                             pastix_bcsc_t  *bcsc )
{
  pastix_int_t  itercol, itercblk;
    pastix_int_t  cblknbr  = solvmtx->cblknbr;
    pastix_int_t  eltnbr   = spm->gNexp;
    pastix_int_t *col2cblk = NULL;

    bcsc->mtxtype = spm->mtxtype;
    bcsc->flttype = spm->flttype;
    bcsc->gN      = spm->gN;
    bcsc->n       = spm->n;

    assert( spm->loc2glob != NULL );

    /*
     * Initialize the col2cblk array. col2cblk[i] contains the cblk index of the
     * i-th column. col2cblk[i] = -1 if not local.
     */
    {
        SolverCblk *cblk = solvmtx->cblktab;

        MALLOC_INTERN( col2cblk, eltnbr, pastix_int_t );
        for (itercol=0; itercol<eltnbr; itercol++)
        {
            col2cblk[itercol] = -1;
        }

        for (itercblk=0; itercblk<cblknbr; itercblk++, cblk++)
        {
            if( cblk->cblktype & (CBLK_FANIN|CBLK_RECV) ){
                continue;
            }
            for (itercol  = cblk->fcolnum;
                 itercol <= cblk->lcolnum;
                 itercol++ )
            {
                col2cblk[itercol] = itercblk;
            }
        }
    }

    /*
     * Fill in the lower triangular part of the blocked csc with values and
     * rows. The upper triangular part is done later if required through LU
     * factorization.
     */
    switch( spm->flttype ) {
    case SpmFloat:
        bcsc_sinit_dist( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmDouble:
        bcsc_dinit_dist( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex32:
        bcsc_cinit_dist( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmComplex64:
        bcsc_zinit_dist( spm, ord, solvmtx, col2cblk, initAt, bcsc );
        break;
    case SpmPattern:
    default:
        fprintf(stderr, "bcsc_init_dist: Error unknown floating type for input spm\n");
    }

    memFree_null(col2cblk);
}

/**
 *******************************************************************************
 *
 * @brief Initialize the block csc matrix.
 *
 * The block csc matrix is used to initialize the factorized matrix, and to
 * perform the matvec operations in refinement.
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
 * @param[in] initAt
 *          A flag to enable/disable the initialization of A'
 *
 * @param[inout] bcsc
 *          On entry, the pointer to an allocated bcsc.
 *          On exit, the bcsc stores the input spm with the permutation applied
 *          and grouped accordingly to the distribution described in solvmtx.
 *
 *******************************************************************************
 *
 * @return The time spent to initialize the bcsc structure.
 *
 *******************************************************************************/
double
bcscInit( const spmatrix_t     *spm,
          const pastix_order_t *ord,
          const SolverMatrix   *solvmtx,
                pastix_int_t    initAt,
                pastix_bcsc_t  *bcsc )
{
    assert( ord->baseval == 0 );
    assert( ord->vertnbr == spm->n );

    double time = 0.;
    clockStart(time);

    if ( spm->loc2glob == NULL ) {
        bcsc_init_centralized( spm, ord, solvmtx, initAt, bcsc );
    }
    else {
        fprintf(stderr, "bcscInit: Distributed SPM not yet supported");
        bcsc_init_dist( spm, ord, solvmtx, initAt, bcsc );
    }

    clockStop(time);
    return time;
}

/**
 *******************************************************************************
 *
 * @brief Cleanup the block csc structure but do not free the bcsc pointer.
 *
 *******************************************************************************
 *
 * @param[inout] bcsc
 *          The block csc matrix to free.
 *
 *******************************************************************************/
void
bcscExit( pastix_bcsc_t *bcsc )
{
    bcsc_cblk_t *cblk;
    pastix_int_t i;

    if ( bcsc->cscftab == NULL ) {
        return;
    }

    for (i=0, cblk=bcsc->cscftab; i < bcsc->cscfnbr; i++, cblk++ ) {
        memFree_null( cblk->coltab );
    }

    memFree_null( bcsc->cscftab );
    memFree_null( bcsc->rowtab );

    if ( (bcsc->Uvalues != NULL) &&
         (bcsc->Uvalues != bcsc->Lvalues) ) {
        memFree_null( bcsc->Uvalues );
    }

    memFree_null( bcsc->Lvalues );
}
