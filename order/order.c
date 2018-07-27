/**
 *
 * @file order.c
 *
 * PaStiX order structure routines
 *
 * @copyright 2004-2018 Bordeaux INP, CNRS (LaBRI UMR 5800), Inria,
 *                      Univ. Bordeaux. All rights reserved.
 *
 * @version 6.0.1
 * @author Francois Pellegrini
 * @author Mathieu Faverge
 * @date 2018-07-16
 *
 **/
#include "common.h"
#include "pastix/order.h"
#include <string.h>

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Allocate the order structure.
 *
 * The base value is set to 0 by default.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call pastixOrderExit() to release the memory first if required to
 *          prevent memory leak.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *          If vertnbr == 0, permtab and peritab are not allocated.
 *
 * @param[in] cblknbr
 *          The number of supernodes. The internal rangtab array is of size
 *          cblknbr+1, and treetab of size cblknbr.
 *          If cblknbr == 0, rangtab and treetab are not allocated.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit,
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect,
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastixOrderAlloc( pastix_order_t * const ordeptr,
                  pastix_int_t           vertnbr,
                  pastix_int_t           cblknbr )
{
    /* Parameter checks */
    if ( ordeptr == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( vertnbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( cblknbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    memset(ordeptr, 0, sizeof(pastix_order_t));

    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;
#if defined(PASTIX_SUPERNODE_STATS)
    ordeptr->sndenbr = cblknbr;
    ordeptr->sndetab = NULL;
#endif

    if (vertnbr != 0) {
        MALLOC_INTERN(ordeptr->permtab, vertnbr, pastix_int_t);
        MALLOC_INTERN(ordeptr->peritab, vertnbr, pastix_int_t);
    }

    if (cblknbr != 0) {
        MALLOC_INTERN(ordeptr->rangtab, cblknbr+1, pastix_int_t);
        MALLOC_INTERN(ordeptr->treetab, cblknbr,   pastix_int_t);
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Initialize the order structure with the given values.
 *
 * The base value is set to 0 by default. This is useful to give a personal
 * ordering to the pastix_task_order() function.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure is set to 0 and then initialize.
 *          Need to call orderExit to release the memory first if required to
 *          prevent memory leak.
 *
 * @param[in] baseval
 *          The base value used in the given arrays. Usually 0 for C, 1 for
 *          Fortran. Must be >= 0.
 *
 * @param[in] vertnbr
 *          The number of nodes, this is the size of the internal permtab and
 *          peritab arrays.
 *
 * @param[in] cblknbr
 *          The number of supernodes. The internal rangtab and treetab arrays
 *          are of size cblknbr+1.
 *
 * @param[in] permtab
 *          The permutation array which must be of size vertnbr, and based on
 *          baseval value.
 *          If NULL, the permtab field is not initialized.
 *
 * @param[in] peritab
 *          The inverse permutation array which must be of size vertnbr, and
 *          based on baseval value.
 *          If NULL, the peritab field is not initialized.
 *
 * @param[in] rangtab
 *          The rangtab array that describes the supernodes in the graph. This
 *          array must be of size cblknbr+1, and based on baseval value.
 *          If NULL, the rangtab field is not initialized.
 *
 * @param[in] treetab
 *          The treetab array that describes the elimination tree connecting the
 *          supernodes. This array must be defined as follow:
 *              - of size cblknbr;
 *              - based on baseval value;
 *              - each treetab[i] > i, unless i is a root and treetab[i] == -1
 *              - all roots of the tree must have -1 as father
 *          If NULL, the treetab field is not initialized.
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 * @retval PASTIX_ERR_OUTOFMEMORY if one allocation failed.
 *
 *******************************************************************************/
int
pastixOrderInit( pastix_order_t * const ordeptr,
                 pastix_int_t           baseval,
                 pastix_int_t           vertnbr,
                 pastix_int_t           cblknbr,
                 pastix_int_t   * const permtab,
                 pastix_int_t   * const peritab,
                 pastix_int_t   * const rangtab,
                 pastix_int_t   * const treetab )
{
    /* Parameter checks */
    if ( ordeptr == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( vertnbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( cblknbr < 0 ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    memset(ordeptr, 0, sizeof(pastix_order_t));

    ordeptr->baseval = baseval;
    ordeptr->vertnbr = vertnbr;
    ordeptr->cblknbr = cblknbr;
#if defined(PASTIX_SUPERNODE_STATS)
    ordeptr->sndenbr = cblknbr;
    ordeptr->sndetab = NULL;
#endif

    if ( permtab ) {
        ordeptr->permtab = permtab;
    }
    if ( peritab ) {
        ordeptr->peritab = peritab;
    }
    if ( rangtab ) {
        ordeptr->rangtab = rangtab;
#if defined(PASTIX_SUPERNODE_STATS)
        ordeptr->sndetab = malloc( (ordeptr->sndenbr+1) * sizeof(pastix_int_t) );
        memcpy( ordeptr->sndetab, ordeptr->rangtab, (ordeptr->sndenbr+1) * sizeof(pastix_int_t) );
#endif
    }
    if ( treetab ) {
        ordeptr->treetab = treetab;
    }

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief Free the arrays initialized in the order structure.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The data structure to clean. All arrays of the structure are freed
 *          and the structure is set to 0.
 *
 *******************************************************************************/
void
pastixOrderExit( pastix_order_t * const ordeptr )
{
    /* Parameter check */
    if ( ordeptr == NULL ) {
        return;
    }

    if (ordeptr->permtab != NULL) {
        memFree_null (ordeptr->permtab);
    }
    if (ordeptr->peritab != NULL) {
        memFree_null (ordeptr->peritab);
    }
    if (ordeptr->rangtab != NULL) {
        memFree_null (ordeptr->rangtab);
    }
    if (ordeptr->treetab != NULL) {
        memFree_null (ordeptr->treetab);
    }
#if defined(PASTIX_SUPERNODE_STATS)
    if (ordeptr->sndetab != NULL) {
        memFree_null (ordeptr->sndetab);
    }
#endif
    memset(ordeptr, 0, sizeof(pastix_order_t) );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine sets the base of the given ordering structure to the
 * given base value.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering to rebase.
 *
 * @param[in] baseval
 *          The base value to be used (needs to be 0 or 1).
 *
 *******************************************************************************/
void
pastixOrderBase( pastix_order_t * const ordeptr,
                 pastix_int_t           baseval )
{
    pastix_int_t baseadj;                    /* Base adjust */
    pastix_int_t cblknum;
    pastix_int_t vertnum;

    /* Parameter checks */
    if ( ordeptr == NULL ) {
        errorPrint("pastixOrderBase: ordeptr pointer is NULL");
        return;
    }
    if ( (baseval != 0) &&
         (baseval != 1) )
    {
        errorPrint("pastixOrderBase: baseval is incorrect, must be 0 or 1");
        return;
    }

    baseadj = baseval - ordeptr->baseval; /* Set base adjust     */
    if (baseadj == 0)                     /* If base already set */
	return;

    if (ordeptr->permtab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++) {
	    ordeptr->permtab[vertnum] += baseadj;
        }
    }
    if (ordeptr->peritab != NULL) {
	for (vertnum = 0; vertnum < ordeptr->vertnbr; vertnum ++) {
	    ordeptr->peritab[vertnum] += baseadj;
        }
    }

    if (ordeptr->rangtab != NULL) {
        for (cblknum = 0; cblknum <= ordeptr->cblknbr; cblknum ++) {
            ordeptr->rangtab[cblknum] += baseadj;
        }
    }
    if (ordeptr->treetab != NULL) {
        for (cblknum = 0; cblknum < ordeptr->cblknbr; cblknum ++) {
            ordeptr->treetab[cblknum] += baseadj;
        }
    }
#if defined(PASTIX_SUPERNODE_STATS)
    if (ordeptr->sndetab != NULL) {
        pastix_int_t sndenum;
        for (sndenum = 0; sndenum <= ordeptr->sndenbr; sndenum ++) {
            ordeptr->sndetab[sndenum] += baseadj;
        }
    }
#endif

    ordeptr->baseval = baseval;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine expand the permutation arrays and the rangtab.
 *
 *******************************************************************************
 *
 * @param[inout] ordeptr
 *          The ordering to expand.
 *
 * @param[in] spm
 *          The sparse matrix structure providing dof information.
 *
 *******************************************************************************/
void pastixOrderExpand( pastix_order_t * const ordeptr,
                        spmatrix_t     * const spm )
{
    pastix_int_t *peritab;
    pastix_int_t  i, j, n, vertnbr;
    pastix_int_t  begin, end;
    pastix_int_t *iter_peri;
    pastix_int_t *rangtab;
    pastix_int_t *dofs;
    pastix_int_t  sum, tmp;

    spmBase( spm, 0 );

    n       = spm->nexp;
    vertnbr = ordeptr->vertnbr;
    MALLOC_INTERN( peritab, vertnbr, pastix_int_t );

    memcpy( peritab, ordeptr->peritab, vertnbr * sizeof(pastix_int_t) );

    ordeptr->peritab = (pastix_int_t*)memRealloc( ordeptr->peritab, n * sizeof(pastix_int_t) );
    ordeptr->permtab = (pastix_int_t*)memRealloc( ordeptr->permtab, n * sizeof(pastix_int_t) );

    iter_peri = ordeptr->peritab;
    /*
     * Initialise permutation tab and its inverse with dofs
     * and previous inverse permutation tab
     */
    sum = 0;
    for (i = 0; i < vertnbr; ++i)
    {
        if ( spm->dof <= 0 ) {
            begin = spm->dofs[ peritab[i] ];
            end   = spm->dofs[ peritab[i] + 1 ];
        }
        else {
            begin = peritab[i] * spm->dof;
            end   = begin + spm->dof;
        }
        sum += end - begin;
        for (j = begin; j < end; ++j, ++iter_peri)
        {
            *iter_peri = j;
            ordeptr->permtab[ j ] = iter_peri - ordeptr->peritab;/* perm[peri[i]] = i*/
        }
    }
    ordeptr->vertnbr = n;

    rangtab = ordeptr->rangtab;
    dofs    = spm->dofs;
    tmp     = 0;
    /*
     * Use previous value of rangtab to compute new ones
     * Use formula :
     * rexp[i] = rexp[i - 1] + sum( r[i], r[i+1] - 1, dofs[peri[j]+1] - dofs[peri[j]] )
     */
    for (i = 1; i <= ordeptr->cblknbr; ++i)
    {
        if (spm->dof > 0) {
            rangtab[i] *= spm->dof ;
        }
        else {
            begin = tmp;        /* r[i - 1] : Updated on last step  */
            end   = rangtab[i]; /* r[i] :     Not yet updated */
            sum   = 0;
            for (j = begin; j < end; ++j)
            {
                sum += dofs[peritab[j] + 1] - dofs[peritab[j]];
            }
            tmp = rangtab[i];
            rangtab[i] = rangtab[i-1] + sum;
        }
    }

    memFree_null( peritab );
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine copy a given ordering in a new one.
 *
 * This function copies an order structure into another one. If all subpointers
 * are NULL, then they are all allocated and conatins the original ordesrc
 * values on exit. If one or more array pointers are not NULL, then, only those
 * are copied to the ordedst structure.
 *
 *******************************************************************************
 *
 * @param[inout] ordedst
 *          The destination ordering
 *
 * @param[in] ordesrc
 *          The source ordering
 *
 *******************************************************************************
 *
 * @retval PASTIX_SUCCESS on successful exit
 * @retval PASTIX_ERR_BADPARAMETER if one parameter is incorrect.
 *
 *******************************************************************************/
int
pastixOrderCopy( pastix_order_t       * const ordedst,
                 const pastix_order_t * const ordesrc )
{
    /* Parameter checks */
    if ( ordedst == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ordesrc == NULL ) {
        return PASTIX_ERR_BADPARAMETER;
    }
    if ( ordesrc == ordedst ) {
        return PASTIX_ERR_BADPARAMETER;
    }

    ordedst->baseval = ordesrc->baseval;
    ordedst->vertnbr = ordesrc->vertnbr;
    ordedst->cblknbr = ordesrc->cblknbr;
#if defined(PASTIX_SUPERNODE_STATS)
    ordedst->sndenbr = ordesrc->sndenbr;
#endif
    if ( (ordedst->permtab == NULL) &&
         (ordedst->peritab == NULL) &&
         (ordedst->rangtab == NULL) &&
         (ordedst->treetab == NULL) )
    {
        pastixOrderAlloc( ordedst, ordedst->vertnbr, ordedst->cblknbr );
    }

    if ( (ordesrc->permtab != NULL) && (ordedst->permtab != NULL) )
    {
        memcpy( ordedst->permtab, ordesrc->permtab, ordesrc->vertnbr * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->peritab != NULL) && (ordedst->peritab != NULL) )
    {
        memcpy( ordedst->peritab, ordesrc->peritab, ordesrc->vertnbr * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->rangtab != NULL) && (ordedst->rangtab != NULL) )
    {
        memcpy( ordedst->rangtab, ordesrc->rangtab, (ordesrc->cblknbr+1) * sizeof(pastix_int_t) );
    }

    if ( (ordesrc->treetab != NULL) && (ordedst->treetab != NULL) )
    {
        memcpy( ordedst->treetab, ordesrc->treetab, ordesrc->cblknbr * sizeof(pastix_int_t) );
    }

#if defined(PASTIX_SUPERNODE_STATS)
    if ( (ordesrc->sndetab != NULL) && (ordedst->sndetab != NULL) )
    {
        memcpy( ordedst->sndetab, ordesrc->sndetab, (ordesrc->sndenbr+1) * sizeof(pastix_int_t) );
    }
#endif
    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_order
 *
 * @brief This routine returns the pointer to the internal order structure to
 * access permutation information.
 *
 * @warning The data returned by the routines must not be freed or modified by
 * the user, and are valid as long as no operation on the ordering is performed
 * (pastix_subtask_order(), pastix_subtask_reordering(), or pastixFinalize()).
 *
 *******************************************************************************
 *
 * @param[in] pastix_data
 *          The pastix_data structure of the problem
 *
 *******************************************************************************
 *
 * @return The pointer to the internal ordering structure with permutation
 *         information, or NULL if pastix_data is invalid.
 *
 *******************************************************************************/
const pastix_order_t *
pastixOrderGet( const pastix_data_t * const pastix_data )
{
    /* Parameter checks */
    if ( pastix_data == NULL ) {
        return NULL;
    }
    return pastix_data->ordemesh;
}
