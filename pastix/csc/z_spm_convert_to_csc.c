/**
 *
 * @file z_spm_convert_to_csc.c
 *
 *  PaStiX csc routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Mathieu Faverge
 * @author Theophile Terraz
 * @date 2015-01-01
 *
 * @precisions normal z -> c d s p
 **/
#include "common.h"
#include "csc.h"

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmConvertIJV2CSC - convert a matrix in IJV format to a matrix in CSC
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The ijv matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertIJV2CSC( pastix_csc_t *spm )
{
#if !defined(PRECISION_p)
    pastix_complex64_t *navals = NULL;
    pastix_complex64_t *oavals = NULL;
#endif
    pastix_int_t       *spmptx, *otmp;
    pastix_int_t i, j, tmp, baseval, total;
    pastix_csc_t oldspm;

    /* Backup the input */
    memcpy( &oldspm, spm, sizeof(pastix_csc_t) );

    /*
     * Check the baseval, we consider that arrays are sorted by columns or rows
     */
    baseval = pastix_imin( *(oldspm.colptr), *(oldspm.rows) );

    /* Compute the new colptr */
    spm->colptr = (pastix_int_t *) calloc(spm->n+1,sizeof(pastix_int_t));

    /* Compute the number of edges per row */
    spmptx = spm->colptr - baseval;
    otmp   = oldspm.colptr;
    for (i=0; i<spm->nnz; i++, otmp++)
    {
        spmptx[ *otmp ] ++;
    }

    /* Compute the indexes in C numbering for the following sort */
    total = 0;
    spmptx = spm->colptr;
    for (i=0; i<(spm->n+1); i++, spmptx++)
    {
        tmp = *spmptx;
        *spmptx = total;
        total += tmp;
    }
    assert( total == spm->nnz );

    /* Sort the rows and avals arrays by column */
    spm->rows  = malloc(spm->nnz * sizeof(pastix_int_t));

#if defined(PRECISION_p)
    spm->avals = NULL;
#else
    spm->avals = malloc(spm->nnz * sizeof(pastix_complex64_t));
    navals = (pastix_complex64_t*)(spm->avals);
    oavals = (pastix_complex64_t*)(oldspm.avals);
#endif

    for (j=0; j<spm->nnz; j++)
    {
        i = oldspm.colptr[j] - baseval;

        spm->rows[ spm->colptr[i] ] = oldspm.rows[j];

#if !defined(PRECISION_p)
        navals[ spm->colptr[i] ] = oavals[j];
#endif
        (spm->colptr[i])++;

        assert( spm->colptr[i] <= spm->colptr[i+1] );
    }

    /* Rebuild the colptr with the correct baseval */
    tmp = spm->colptr[0];
    spm->colptr[0] = baseval;

    spmptx = spm->colptr + 1;
    for (i=1; i<(spm->n+1); i++, spmptx++)
    {
        total = *spmptx;
        *spmptx = tmp + baseval;
        tmp = total;
    }
    assert( spm->colptr[ spm->n ] == (spm->nnz+baseval) );

    free( oldspm.colptr );
    free( oldspm.rows );

    if (oldspm.avals != NULL)
        free( oldspm.avals );

    spm->fmttype = PastixCSC;

    return PASTIX_SUCCESS;
}

/**
 *******************************************************************************
 *
 * @ingroup pastix_csc
 *
 * z_spmConvertCSR2CSC - convert a matrix in CSR format to a matrix in CSC
 * format.
 *
 *******************************************************************************
 *
 * @param[in,out] spm
 *          The csr matrix at enter,
 *          the csc matrix at exit.
 *
 *******************************************************************************
 *
 * @return
 *      \retval PASTIX_SUCCESS
 *
 *******************************************************************************/
int
z_spmConvertCSR2CSC( pastix_csc_t *spm )
{
    assert( spm->loc2glob == NULL );
    assert( spm->fmttype == PastixCSR );

    spm->fmttype = PastixCSC;

    switch( spm->mtxtype ) {
#if defined(PRECISION_z) || defined(PRECISION_c)
    case PastixHermitian:
    {
        /* Similar to PastixSymmetric case with conjugate of the values */
        pastix_complex64_t *valptr = spm->avals;
        pastix_int_t i;

        for(i=0; i<spm->nnz; i++, valptr++){
            *valptr = conj( *valptr );
        }
    }
#endif
    case PastixSymmetric:
    {
        pastix_int_t *tmp;

        /* Just need to swap the pointers */
        tmp         = spm->rows;
        spm->rows   = spm->colptr;
        spm->colptr = tmp;

        return PASTIX_SUCCESS;
    }
    break;

    case PastixGeneral:
    default:
    {
        pastix_int_t       *row_csc;
        pastix_int_t       *col_csc;
#if !defined(PRECISION_p)
        pastix_complex64_t *val_csc;
        pastix_complex64_t *valptr = (pastix_complex64_t*)(spm->avals);
#endif
        pastix_int_t j, k, col, row, nnz, baseval;

        baseval = pastix_imin( *(spm->colptr), *(spm->rows) );
        nnz = spm->nnz;

        row_csc = malloc(nnz * sizeof(pastix_int_t));
        col_csc = calloc(spm->n+1,sizeof(pastix_int_t));

        assert( row_csc );
        assert( col_csc );

#if !defined(PRECISION_p)
        val_csc = malloc(nnz*sizeof(pastix_complex64_t));
        assert( val_csc );
#endif

        /* Count the number of elements per column */
        for (j=0; j<nnz; j++) {
            col = spm->colptr[j] - baseval;
            col_csc[ col+1 ] ++;
        }

        /* Compute the index of each column */
        col_csc[0] = 0;
        for (j=0; j<spm->n; j++){
            col_csc[j+1] += col_csc[j];
        }

        assert( (col_csc[spm->gN]-baseval) == nnz );

        for (row=0; row<spm->n; row++) {
            pastix_int_t fcol = spm->rows[row  ] - baseval;
            pastix_int_t lcol = spm->rows[row+1] - baseval;

            for (k=fcol; k<lcol; k++) {
                col = spm->colptr[k] - baseval;
                j = col_csc[col];
                row_csc[j] = row + baseval;

#if !defined(PRECISION_p)
                val_csc[j] = valptr[k];
#endif
                col_csc[col] ++;
            }
        }

        /* Restore the colptr indexes */
        {
            pastix_int_t tmp, tmp2;

            tmp = col_csc[0];
            col_csc[0] = baseval;
            for (j=0; j<spm->n; j++) {
                tmp2 = col_csc[j+1];
                col_csc[j+1] = tmp;
                tmp = tmp2;
            }
        }

        memFree_null(spm->colptr);
        memFree_null(spm->rows);
        spm->colptr = col_csc;
        spm->rows   = row_csc;
#if !defined(PRECISION_p)
        memFree_null(spm->avals);
        spm->avals = val_csc;
#else
        spm->avals = NULL;
#endif
    }
    }

    return PASTIX_SUCCESS;
}
