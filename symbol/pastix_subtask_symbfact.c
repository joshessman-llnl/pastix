/**
 *
 * @file pastix_subtask_symbfact.c
 *
 * PaStiX symbolic factorizations task.
 * Contains wrappers to the symbolic factorization step.
 * Affetcted by the compilation time options:
 *    - PASTIX_SYMBOL_DUMP_SYMBMTX: Dump the symbol matrix in a postscript file.
 *    - COMPACT_SMX: Optimization for solve step (TODO: check if not obsolete)
 *    - FORGET_PARTITION: Force to forget the precomputed partition
 *
 * @copyright 2015-2020 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.3
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2019-11-12
 *
 **/
#include "common.h"
#include "spm.h"
#include "graph.h"
#include "pastix/order.h"
#include "symbol.h"
#if defined( PASTIX_DISTRIBUTED )
#include "csc_utils.h"
#include "cscd_utils_intern.h"
#endif /* defined(PASTIX_DISTRIBUTED) */

#define searchInList(n, list, list_size, index){  \
    pastix_int_t min = 0;                                  \
    pastix_int_t max = list_size-1;                        \
    if  (list_size == 0) {                        \
      index = -1;                                 \
    }                                             \
    else {                                        \
      while (1) {                                 \
        pastix_int_t cursor = (max +min)/2;                \
        if ((list)[cursor] == n)                  \
          {                                       \
            index = cursor;                       \
            break;                                \
          }                                       \
        else {                                    \
          if (max <= min) {                       \
            index = -1;                           \
            break;                                \
          }                                       \
          if ((list)[cursor] > n)                 \
            {                                     \
              max = cursor-1;                     \
            }                                     \
          else {                                  \
            if ((list)[cursor] < n)               \
              {                                   \
                min = cursor+1;                   \
              }                                   \
          }                                       \
        }                                         \
      }                                           \
    }                                             \
  }

/*
 *   Function: cscd2csc_int
 *
 *   Transform a cscd to a csc.
 *   colptr2, row2, avals2, rhs2, perm2, invp2 are allocated here.
 *
 *   Parameters:
 *     lN          - number of local column.
 *     lcolptr     - starting index of each local column in row and avals.
 *     lrow        _ row number of each local element.
 *     lavals      - values of each local element.
 *     lrhs        - local part of the right hand side.
 *     lperm       - local part of the permutation tabular.
 *     linvp       - Means nothing, to suppress.
 *     gN          - global number of columns (output).
 *     gcolptr     - starting index of each column in row2 and avals2 (output).
 *     grow        - row number of each element (output).
 *     gavals      - values of each element (output).
 *     grhs        - global right hand side (output).
 *     gperm       - global permutation tabular (output).
 *     ginvp       - global reverse permutation tabular (output).
 *     loc2glob    - global number of each local column.
 *     pastix_comm - PaStiX MPI communicator.
 *     intern_flag - Decide if malloc will use internal or external macros.
 *
 */
void  cscd2csc_int(pastix_int_t  lN, pastix_int_t *  lcolptr, pastix_int_t * lrow, double * lavals,
                   double * lrhs, pastix_int_t * lperm, pastix_int_t * linvp,
                   pastix_int_t *gN, pastix_int_t ** gcolptr, pastix_int_t **grow, double **gavals,
                   double **grhs, pastix_int_t **gperm, pastix_int_t **ginvp,
                   pastix_int_t *loc2glob, PASTIX_Comm pastix_comm, pastix_int_t ndof, int intern_flag)
{
  pastix_int_t      i,j;
  int      rank;
  int      commSize;
  pastix_int_t      index;
  pastix_int_t    * AllLocalN   = NULL;
  pastix_int_t   ** AllColptr   = NULL;
  pastix_int_t   ** AllRow      = NULL;
  pastix_int_t   ** AllLoc2glob = NULL;
  double ** AllAvals    = NULL;
  double ** AllRhs      = NULL;
  pastix_int_t   ** AllPerm     = NULL;
  /*   pastix_int_t   ** AllInvp; */
  int      proc;
  (void)pastix_comm; (void)linvp;

  MPI_Comm_rank(pastix_comm,&rank);
  MPI_Comm_size(pastix_comm,&commSize);

  MALLOC_INTERN(AllLocalN,   commSize, pastix_int_t );
  MALLOC_INTERN(AllColptr,   commSize, pastix_int_t*);
  MALLOC_INTERN(AllRow,      commSize, pastix_int_t*);
  MALLOC_INTERN(AllLoc2glob, commSize, pastix_int_t*);

  if (gavals != NULL)
    MALLOC_INTERN(AllAvals, commSize, double*);

  if (grhs != NULL)
    MALLOC_INTERN(AllRhs, commSize, double*);

  if (gperm != NULL)
    MALLOC_INTERN(AllPerm, commSize, pastix_int_t*);

  for (i = 0; i < commSize; i++)
    {
      if (i == rank)
        {
          AllLocalN[i]   = lN;
          AllColptr[i]   = lcolptr;
          AllRow[i]      = lrow;
          if (gavals != NULL)
            AllAvals[i]    = lavals;
          if (grhs != NULL)
            AllRhs[i]      = lrhs;
          if (gperm != NULL)
            AllPerm[i]     = lperm;
          AllLoc2glob[i] = loc2glob;
        }

      MPI_Bcast(&AllLocalN[i] , 1, PASTIX_MPI_INT, i, pastix_comm);
      if (rank != i)
        {
          MALLOC_INTERN(AllColptr[i],   AllLocalN[i]+1, pastix_int_t);
          MALLOC_INTERN(AllLoc2glob[i], AllLocalN[i],   pastix_int_t);
          if (grhs != NULL)
            MALLOC_INTERN(AllRhs[i], ndof*AllLocalN[i], double);

          if (gperm != NULL)
            MALLOC_INTERN(AllPerm[i], AllLocalN[i], pastix_int_t);
        }

      MPI_Bcast(AllColptr[i], AllLocalN[i]+1, PASTIX_MPI_INT  , i, pastix_comm);
      if (rank != i)
        {
          MALLOC_INTERN(AllRow[i], AllColptr[i][AllLocalN[i]]-1, pastix_int_t);
          if (gavals != NULL)
            MALLOC_INTERN(AllAvals[i], ndof*ndof*(AllColptr[i][AllLocalN[i]]-1), double);
        }

      MPI_Bcast(AllRow[i], AllColptr[i][AllLocalN[i]]-1,
                PASTIX_MPI_INT, i, pastix_comm);
      MPI_Bcast(AllLoc2glob[i], AllLocalN[i], PASTIX_MPI_INT, i, pastix_comm);
      if (gperm != NULL) {
        MPI_Bcast(AllPerm[i], AllLocalN[i], PASTIX_MPI_DOUBLE, i, pastix_comm);
      }

      if (grhs != NULL) {
        MPI_Bcast(AllRhs[i], ndof*AllLocalN[i], PASTIX_MPI_DOUBLE, i, pastix_comm);
      }
      if (gavals != NULL) {
        MPI_Bcast(AllAvals[i], ndof*ndof*(AllColptr[i][AllLocalN[i]]-1),
                  PASTIX_MPI_DOUBLE, i, pastix_comm);
      }
    }


  *gN = 0;
  for (i = 0; i < commSize; i++)
    {
      *gN += AllLocalN[i];
    }

  MALLOC_INTOREXTERN(*gcolptr, ((*gN)+1), pastix_int_t, intern_flag);
  /* Fill in gcolptr */
  (*gcolptr)[0] = 1;
  for (i = 0; i < (*gN); i++)
    {
      for (proc = 0; proc < commSize; proc++)
        {
            // Check for baseval 0
            if (AllLoc2glob[0][0] == 0) {
                for (j = 0; j < AllLocalN[proc]; j++)
                {
                    AllLoc2glob[proc][j] += 1;
                }
            }
          searchInList(i+1, AllLoc2glob[proc], AllLocalN[proc], index);
          if ( index >= 0 )
            {
              (*gcolptr)[i+1] = (*gcolptr)[i] +
                AllColptr[proc][index+1] -
                AllColptr[proc][index];
              break;
            }

        }
    }
  MALLOC_INTOREXTERN(*grow, (*gcolptr)[(*gN)]-1, pastix_int_t, intern_flag);
  if (gavals != NULL)
    MALLOC_INTOREXTERN(*gavals, ndof*ndof*((*gcolptr)[(*gN)]-1), double, intern_flag);
  if (grhs != NULL)
    MALLOC_INTOREXTERN(*grhs, *gN*ndof, double, intern_flag);
  if (gperm != NULL)
    {
      MALLOC_INTOREXTERN(*gperm, *gN, pastix_int_t, intern_flag);
      MALLOC_INTOREXTERN(*ginvp, *gN, pastix_int_t, intern_flag);
    }

  /* Fill-in grow, gavals, grhs, gperm and ginvp*/
  for (proc = 0; proc < commSize; proc++)
    {
      for (i = 0; i < AllLocalN[proc]; i++)
        {
          memcpy(&(*grow)[(*gcolptr)[AllLoc2glob[proc][i]-1]-1],
                 &AllRow[proc][AllColptr[proc][i]-1],
                 (AllColptr[proc][i+1] - AllColptr[proc][i])*sizeof(pastix_int_t));
          if (gavals != NULL)
            memcpy(&(*gavals)[(*gcolptr)[AllLoc2glob[proc][i]-1]-1],
                   &AllAvals[proc][AllColptr[proc][i]-1],
                   ndof*ndof*(AllColptr[proc][i+1] - AllColptr[proc][i])*sizeof(double));
          if (grhs != NULL)
            for (j = 0; j < ndof; j++)
              (*grhs)[ndof*(AllLoc2glob[proc][i]-1)+j]  = AllRhs[proc][ndof*i+j];
          if (gperm != NULL)
            {
              (*gperm)[AllLoc2glob[proc][i]-1] = AllPerm[proc][i];
              /* (*ginvp)[AllLoc2glob[proc][i]-1] = AllInvp[proc][i]; */

            }
        }
    }
  if (gperm != NULL)
    {
      for (i = 0; i < *gN; i++)
        (*ginvp)[(*gperm)[i]] = i;
    }
  for (i = 0; i < commSize; i++)
    {
      if (rank != i)
        {
          memFree_null(AllColptr[i]);
          memFree_null(AllRow[i]);
          if (gavals != NULL)
            memFree_null(AllAvals[i]);
          if (grhs != NULL)
            memFree_null(AllRhs[i]);
          memFree_null(AllLoc2glob[i]);
          if (gperm != NULL)
            {
              memFree_null(AllPerm[i]);
              /*            memFree_null(AllInvp[i]); */
            }
        }
    }
  memFree_null(AllLocalN);
  memFree_null(AllColptr);
  memFree_null(AllRow);
  if (gavals != NULL)
    memFree_null(AllAvals);
  if (grhs != NULL)
    memFree_null(AllRhs);
  memFree_null(AllLoc2glob);
  if (gperm != NULL)
    {
      memFree_null(AllPerm);
      /*       memFree_null(AllInvp); */
    }

}

/**
 *******************************************************************************
 *
 * @ingroup pastix_analyze
 *
 * @brief Computes the symbolic factorization step.
 *
 * Computes the symbolic matrix structure and if required the amalgamated
 * supernode partition.
 *
 * The function is a *centralized* algorithm to generate the symbol matrix
 * structure associated to the problem. It takes as input the ordemesh structure
 * (permutation array, inverse permutation array, and optionnal supernodes
 * array) and returns the modified ordemesh structure if changed, and the
 * symbolic structure.
 *  - If a direct factorization is performed, the structure is generated with
 * pastixSymbolFaxDirect() thanks to the information provided by the ordering
 * steps (permutation, partition, and elimination tree).
 * - If an ILU(k) factorization is performed, pastixSymbolFaxILUk() is used to
 * generate the symbol matrix structure. It requires an intermediate step to
 * generate the csr graph of L with incomplete factorization.
 *
 * Both algorithms are working with a centralized version of the graph and are
 * replicated on every nodes. If a distributed graph has been used, it is
 * gathered on each node to compute the symbol matrix.
 *
 * This routine is affected by the following parameters:
 *   IPARM_VERBOSE, IPARM_INCOMPLETE, IPARM_LEVEL_OF_FILL,
 *   IPARM_IO_STRATEGY, IPARM_FLOAT, IPARM_FACTORIZATION
 *
 * On exit, the following parameters are set:
 *   IPARM_NNZEROS, DPARM_FACT_THFLOPS, DPARM_FACT_RLFLOPS
 *
 *******************************************************************************
 *
 * @param[inout] pastix_data
 *          The pastix_data structure that describes the solver instance.
 *          On exit, the field symbmtx is initialized with the symbol matrix,
 *          and the field ordemesh is updated if the supernode partition is
 *          computed.
 *          - IPARM_INCOMPLETE switches the factorization mode from direct to ILU(k).
 *          - IPARM_LEVEL_OF_FILL defines the level of incomplete factorization
 *            if IPARM_INCOMPLETE == 1. If IPARM_LEVEL_OF_FILL < 0, the
 *            full pattern is generated as for direct factorization.
 *          - IPARM_IO_STRATEGY will enable to load/store the result to files.
 *          If set to PastixIOSave, the symbmtx and the generated ordemesh is
 *          dump to file.
 *          If set to PastixIOLoad, the symbmtx (only) is loaded from the files.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 * @retval PASTIX_ERR_INTEGER_TYPE if Scotch integer type is not the
 *         same size as PaStiX ones.
 * @retval PASTIX_ERR_INTERNAL if an error occurs internally to Scotch.
 *
 *******************************************************************************/
int
pastix_subtask_symbfact( pastix_data_t *pastix_data )
{
    pastix_int_t   *iparm;
    double         *dparm;
    pastix_graph_t *graph;
    pastix_order_t *ordemesh;
    int             procnum;
    Clock           timer;

#if defined( PASTIX_DISTRIBUTED )
    pastix_int_t *PTS_perm     = NULL;
    pastix_int_t *PTS_rev_perm = NULL;
    pastix_int_t *tmpperm      = NULL;
    pastix_int_t *tmpperi      = NULL;
    pastix_int_t  gN;
    pastix_int_t  i;
#endif

    /*
     * Check parameters
     */
    if ( pastix_data == NULL ) {
        errorPrint( "pastix_subtask_symbfact: wrong pastix_data parameter" );
        return PASTIX_ERR_BADPARAMETER;
    }
    iparm = pastix_data->iparm;
    dparm = pastix_data->dparm;

    if ( !( pastix_data->steps & STEP_ORDERING ) ) {
        errorPrint( "pastix_subtask_symbfact: pastix_subtask_order() has to be called before "
                    "calling this function" );
        return PASTIX_ERR_BADPARAMETER;
    }

    procnum  = pastix_data->procnum;
    graph    = pastix_data->graph;
    ordemesh = pastix_data->ordemesh;

    if ( graph == NULL ) {
        errorPrint( "pastix_subtask_symbfact: the pastix_data->graph field has not been "
                    "initialized, pastix_subtask_order should be called first" );
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ordemesh == NULL ) {
        errorPrint( "pastix_subtask_symbfact: the pastix_data->ordemesh field has not been "
                    "initialized, pastix_subtask_order should be called first" );
        return PASTIX_ERR_BADPARAMETER;
    }

    clockStart( timer );

    /* Make sure they are both 0-based */
    pastixOrderBase( ordemesh, 0 );
    graphBase( graph, 0 );

    if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
        pastix_print( procnum, 0, OUT_STEP_FAX );
    }

    /* Allocate the symbol matrix structure */
    if ( pastix_data->symbmtx == NULL ) {
        MALLOC_INTERN( pastix_data->symbmtx, 1, symbol_matrix_t );
    }
    else {
        pastixSymbolExit( pastix_data->symbmtx );
    }

    /*Symbol matrix loaded from file */
    if ( iparm[IPARM_IO_STRATEGY] & PastixIOLoad ) {
        FILE *stream = NULL;
        stream       = pastix_fopen( "symbname" );
        if ( stream ) {
            pastixSymbolLoad( pastix_data->symbmtx, stream );
            fclose( stream );
        }
    }
    /* Symbol matrix computed through Fax (Direct or ILU(k)) */
    else {
        pastix_int_t  nfax;
        pastix_int_t *colptrfax;
        pastix_int_t *rowfax;

        /* Check correctness of parameters */
        if ( iparm[IPARM_INCOMPLETE] == 0 ) {
#if defined( COMPACT_SMX )
            if ( procnum == 0 )
                errorPrintW( "COMPACT_SMX only works with incomplete factorization, force ILU(%d) "
                             "factorization.",
                             iparm[IPARM_LEVEL_OF_FILL] );
            iparm[IPARM_INCOMPLETE] = 1;
#endif
        }
        /* End of parameters check */

        /*
         * Fax works with centralized interface, we convert the cscd to csc if required
         */
// #if defined( PASTIX_DISTRIBUTED )
#if defined(PASTIX_WITH_MPI)
        if ( graph->loc2glob != NULL ) {
            cscd2csc_int( graph->n,
                          graph->colptr,
                          graph->rows,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          &nfax,
                          &colptrfax,
                          &rowfax,
                          NULL,
                          NULL,
                          NULL,
                          NULL,
                          graph->loc2glob,
                          pastix_data->pastix_comm,
                          iparm[IPARM_DOF_NBR],
                          1 );
        }
        else
#endif
        {
            nfax      = graph->n;
            colptrfax = graph->colptr;
            rowfax    = graph->rows;
        }

        pastixSymbolInit( graph, ordemesh, pastix_data->symbmtx );

        /*
         * The amalgamate supernodes partition has been found with (PT-)Scotch,
         * we use it to generate the symbol matrix structure.
         * This works only if direct factorization will be performed.
         */
        if ( !iparm[IPARM_INCOMPLETE] || ( iparm[IPARM_LEVEL_OF_FILL] == -1 ) ) {
            if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
                pastix_print( procnum, 0, OUT_FAX_METHOD, "Fax Direct" );
            }

            pastixSymbolFaxDirect( pastix_data->symbmtx, /* Symbol Matrix   */
                                   graph,
                                   ordemesh );
        }
        else {
            if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
                pastix_print( procnum, 0, OUT_FAX_METHOD, "Fax ILU(k)" );
            }

            pastixSymbolFaxILUk( pastix_data->symbmtx, /* Symbol Matrix   */
                                 iparm[IPARM_LEVEL_OF_FILL],
                                 graph,
                                 ordemesh );
        }

        /* Set the beginning of the Schur complement */
        pastix_data->symbmtx->schurfcol =
            nfax - pastix_data->schur_n + pastix_data->symbmtx->baseval;

        if ( graph->loc2glob != NULL ) {
            memFree_null( colptrfax );
            memFree_null( rowfax );
        }
    } /* not PastixIOLoad */

    /* Rebase to 0 */
    pastixSymbolBase( pastix_data->symbmtx, 0 );

    /* Build the browtabs and Realign data structure */
    pastixSymbolBuildRowtab( pastix_data->symbmtx );
    pastixSymbolRealloc( pastix_data->symbmtx );

    if ( ordemesh->selevtx != NULL ) {
        symbol_matrix_t *symbmtx = pastix_data->symbmtx;
        symbol_cblk_t   *cblk    = symbmtx->cblktab;
        int8_t          *selevtx = ordemesh->selevtx;
        pastix_int_t i;

        for(i=0; i<symbmtx->cblknbr; i++, cblk++, selevtx++ ) {
            cblk->selevtx = *selevtx;
        }
    }

#if !defined( NDEBUG )
    if ( pastixOrderCheck( ordemesh ) != 0 ) {
        errorPrint( "pastix_subtask_symbfact: pastixOrderCheck on final ordering after symbolic "
                    "factorization failed !!!" );
        assert( 0 );
    }
    if ( pastixSymbolCheck( pastix_data->symbmtx ) != 0 ) {
        errorPrint( "pastix_subtask_symbfact: symbolCheck on final symbol matrix failed !!!" );
        assert( 0 );
    }
#endif

    /*
     * Save the symbolic factorization
     */
    if ( iparm[IPARM_IO_STRATEGY] & PastixIOSave ) {
        pastix_gendirectories( pastix_data );
        if ( procnum == 0 ) {
            FILE *stream = NULL;
            stream       = pastix_fopenw( pastix_data->dir_global, "symbgen", "w" );
            if ( stream ) {
                pastixSymbolSave( pastix_data->symbmtx, stream );
                fclose( stream );
            }
        }
    }

    /*
     * Dump an eps file of the symbolic factorization
     */
#if defined( PASTIX_SYMBOL_DUMP_SYMBMTX )
    {
        pastix_gendirectories( pastix_data );
        if ( procnum == 0 ) {
            FILE *stream = NULL;
            stream       = pastix_fopenw( pastix_data->dir_global, "symbol.eps", "w" );
            if ( stream ) {
                pastixSymbolDraw( pastix_data->symbmtx, stream );
                fclose( stream );
            }
        }
    }
#endif

    /*
     * Computes statistics and print informations
     */
    iparm[IPARM_NNZEROS] = pastixSymbolGetNNZ( pastix_data->symbmtx );
    pastixSymbolGetFlops( pastix_data->symbmtx,
                          iparm[IPARM_FLOAT],
                          iparm[IPARM_FACTORIZATION],
                          &( dparm[DPARM_FACT_THFLOPS] ),
                          &( dparm[DPARM_FACT_RLFLOPS] ) );

    clockStop( timer );

    /* Warning: the timer will be overwritten by the call in reordering step */
    pastix_data->dparm[DPARM_SYMBFACT_TIME] = clockVal(timer);

    if ( procnum == 0 ) {
        if ( iparm[IPARM_VERBOSE] > PastixVerboseNo )
            pastixSymbolPrintStats( pastix_data->symbmtx );

        if ( iparm[IPARM_VERBOSE] > PastixVerboseNot ) {
            double fillin =
                (double)( iparm[IPARM_NNZEROS] ) / (double)( ( pastix_data->csc )->gnnz );

            pastix_print( procnum, 0, OUT_FAX_SUMMARY,
                          (long)iparm[IPARM_NNZEROS],
                          fillin, clockVal(timer) );
        }
    }

    /* Invalidate following steps, and add order step to the ones performed */
    pastix_data->steps &= ~( STEP_ANALYSE   |
                             STEP_CSC2BCSC  |
                             STEP_BCSC2CTAB |
                             STEP_NUMFACT   |
                             STEP_SOLVE     |
                             STEP_REFINE    );
    pastix_data->steps |= STEP_SYMBFACT;

    return PASTIX_SUCCESS;
}
