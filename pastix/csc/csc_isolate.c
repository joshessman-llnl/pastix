/**
 *
 * @file csc_isolate.h
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 **/
#include "common.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * cscIsolate - This routine isolate a subset of vertices from a given csc, and
 * return a new CSC cleaned from those vertices.
 *
 *******************************************************************************
 *
 * @param[in] n
 *          The number of columns of the original CSC matrix.
 *
 * @param[in] colptr
 *          Array of size n+1
 *          Index of first element of each column in rows array.
 *
 * @param[in] rows
 *          Array of size nnz = colptr[n] - colptr[0].
 *          Rows of each non zero entries.
 *
 * @param[in] loc2glob
 *          Array of size n
 *          Global numbering of each local vertex.
 *
 * @param[in] isolate_n
 *          The number of columns to isolate from the original CSC matrix.
 *
 * @param[in,out] isolate_list
 *          Array of size isolate_n.
 *          List of columns to isolate. On exit, the list is sorted by ascending
 *          indexes.
 *
 * @param[out] new_colptr
 *          Array of size n-isolate_n+1
 *          Index of first element of each column in rows array for the new CSC
 *          matrix.
 *          If new_colptr == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_rows
 *          Array of size new_nnz = (*new_colptr)[n] - (*new_colptr)[0].
 *          Rows of each non zero entries for the new CSC matrix.
 *          If new_rows == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_perm
 *          Array of size n-isolate_n.
 *          Contains permutation generated to isolate the columns at the end of
 *          the CSC.
 *          If new_perm == NULL, nothing is returned, otherwise the pointer to
 *          the allocated sturcture.
 *
 * @param[out] new_invp
 *          Array of size n-isolate_n.
 *          Contains the inverse permutation genereated to isolate the columns
 *          at the end of the CSC.
 *          If new_invp == NULL, nothing is returned, otherwise the pointer to
 *          the allocated structure.
 *
 * @param[out] new_rows
 *          Array of size new_nnz = (*new_colptr)[n] - (*new_colptr)[0].
 *          Rows of each non zero entries for the new CSC matrix.
 *
 * @param[out] loc2glob
 *          Array of size n
 *          Global numbering of each local vertex.
 *
 *******************************************************************************
 *
 * @return
 *          \retval 0 on success.
 *          \retval !0 on failure.
 *
 *******************************************************************************/
/*
 *    perm         - permutation tabular.
 *    invp         - reverse permutation tabular.
*/
int cscIsolate(       pastix_int_t   n,
                const pastix_int_t  *colptr,
                const pastix_int_t  *rows,
                      pastix_int_t   isolate_n,
                      pastix_int_t  *isolate_list,
                      pastix_int_t **new_colptr,
                      pastix_int_t **new_rows,
                      pastix_int_t **new_perm,
                      pastix_int_t **new_invp )
{
    pastix_int_t *tmpcolptr = NULL;
    pastix_int_t *tmprows   = NULL;
    pastix_int_t *tmpperm   = NULL;
    pastix_int_t *tmpinvp   = NULL;
    pastix_int_t  baseval = colptr[0];
    pastix_int_t  nnz = colptr[n] - baseval;
    pastix_int_t  new_n = n - isolate_n;
    pastix_int_t  new_nnz;
    pastix_int_t  i, j, ip, k;
    pastix_int_t  iter_isolate = 0;
    pastix_int_t  iter_non_isolate  = 0;

    if (isolate_n > n) {
        errorPrintW( "Number of columns to isolate greater than the columns in the CSC matrix\n");
        return PASTIX_ERR_BADPARAMETER;
    }

    /* Quick Return */
    if (isolate_n == 0) {
        if (new_colptr != NULL) *new_colptr = (pastix_int_t*)colptr;
        if (new_rows   != NULL) *new_rows   = (pastix_int_t*)rows;
        return PASTIX_SUCCESS;
    }

    if (isolate_n == n) {
        if (new_colptr != NULL) {
            MALLOC_INTERN(*new_colptr, n, pastix_int_t);
            memcpy( *new_colptr, colptr, n*sizeof(pastix_int_t) );
        }
        if (new_rows != NULL) {
            MALLOC_INTERN(*new_rows, nnz, pastix_int_t);
            memcpy( *new_rows, rows, nnz*sizeof(pastix_int_t) );
        }
        return PASTIX_SUCCESS;
    }

    /* Sort the lost of vertices */
    intSort1asc1(isolate_list, isolate_n);

    /* Init invp array */
    MALLOC_INTERN(tmpinvp, n, pastix_int_t);
    for (i = 0; i <n; i++) {
        if ((iter_isolate < isolate_n)          &&
            (i == isolate_list[iter_isolate]-baseval) )
        {
            tmpinvp[new_n+iter_isolate] = i;
            iter_isolate++;
        }
        else
        {
            tmpinvp[iter_non_isolate] = i;
            iter_non_isolate++;
        }
    }

    assert(iter_non_isolate == new_n    );
    assert(iter_isolate     == isolate_n);

    /* Init perm array */
    MALLOC_INTERN(tmpperm, n, pastix_int_t);
    for(i = 0; i < n; i++)
        tmpperm[tmpinvp[i]] = i;

#if defined(PASTIX_DEBUG_CSC)
    for(i = 0; i < n; i++)
    {
        assert(tmpperm[i] < n );
        assert(tmpperm[i] > -1);
    }
#endif

    /* Create the new_colptr array */
    MALLOC_INTERN(tmpcolptr, new_n + 1, pastix_int_t);
    memset(tmpcolptr, 0, (new_n + 1)*sizeof(pastix_int_t));

    tmpcolptr[0] = baseval;
    for (i=0; i<n; i++)
    {
        ip = tmpperm[i];
        if (ip < new_n)
        {
            for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j++)
            {
                /* Count edges in each column of the new graph */
                if (tmpperm[rows[j]-baseval] < new_n)
                {
                    tmpcolptr[ip+1]++;
                }
            }
        }
    }

    for (i = 0; i < new_n; i++)
        tmpcolptr[i+1] += tmpcolptr[i];

    new_nnz = tmpcolptr[new_n] - tmpcolptr[0];
    assert( new_nnz > new_n );

    /* Create the new rows array */
    MALLOC_INTERN(tmprows, new_nnz, pastix_int_t);
    for (i = 0; i <n; i++)
    {
        ip = tmpperm[i];
        if (ip < new_n)
        {
            k = tmpcolptr[ip]-baseval;
            for (j = colptr[i]-baseval; j < colptr[i+1]-baseval; j ++)
            {
                /* Count edges in each column of the new graph */
                if (tmpperm[rows[j]-baseval] < new_n)
                {
                    tmprows[k] = tmpperm[rows[j]-baseval] + baseval;
                    k++;
                }
            }
            assert( k == tmpcolptr[ip+1]-baseval );
        }
    }

    if (new_colptr != NULL) {
        *new_colptr = tmpcolptr;
    } else {
        memFree_null( tmpcolptr );
    }
    if (new_rows != NULL) {
        *new_rows = tmprows;
    } else {
        memFree_null( tmprows );
    }
    if (new_perm != NULL) {
        *new_perm = tmpperm;
    } else {
        memFree_null( tmpperm );
    }
    if (new_invp != NULL) {
        *new_invp = tmpinvp;
    } else {
        memFree_null( tmpinvp );
    }

    return PASTIX_SUCCESS;
}
