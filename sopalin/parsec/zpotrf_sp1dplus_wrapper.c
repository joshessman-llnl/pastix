/*
 * Copyright (c) 2010      The University of Tennessee and The University
 *                         of Tennessee Research Foundation.  All rights
 *                         reserved.
 *
 * @precisions normal z -> s d c
 *
 */
#include <parsec.h>
#include <parsec/data_distribution.h>
#include <parsec/private_mempool.h>
#include <parsec/arena.h>
#include <data_dist/matrix/matrix.h>
#include "common.h"
#include "solver.h"
#include "sopalin_data.h"
#include "parsec/zpotrf_sp1dplus.h"

parsec_handle_t*
dsparse_zpotrf_sp_New( sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    parsec_zpotrf_sp1dplus_handle_t *parsec_zpotrf_sp = NULL;

    parsec_zpotrf_sp = parsec_zpotrf_sp1dplus_new( A, sopalin_data, NULL );

    parsec_zpotrf_sp->_g_p_work = (parsec_memory_pool_t*)malloc(sizeof(parsec_memory_pool_t));
    parsec_private_memory_init( parsec_zpotrf_sp->_g_p_work, sopalin_data->solvmtx->gemmmax * sizeof(pastix_complex64_t) );

    parsec_matrix_add2arena_rect( parsec_zpotrf_sp->arenas[PARSEC_zpotrf_sp1dplus_DEFAULT_ARENA],
                                  parsec_datatype_double_complex_t,
                                  /*sopalin_data->solvmtx->gemmmax*/ 1, 1, 1 );

    return (parsec_handle_t*)parsec_zpotrf_sp;
}

void
dsparse_zpotrf_sp_Destruct( parsec_handle_t *handle )
{
    parsec_zpotrf_sp1dplus_handle_t *parsec_zpotrf_sp = NULL;
    parsec_zpotrf_sp = (parsec_zpotrf_sp1dplus_handle_t *)handle;

    parsec_matrix_del2arena( parsec_zpotrf_sp->arenas[PARSEC_zpotrf_sp1dplus_DEFAULT_ARENA] );

    parsec_private_memory_fini( parsec_zpotrf_sp->_g_p_work );
    free( parsec_zpotrf_sp->_g_p_work );

    parsec_handle_free( handle );
}

int dsparse_zpotrf_sp( parsec_context_t *parsec,
                       sparse_matrix_desc_t *A,
                       sopalin_data_t *sopalin_data )
{
    parsec_handle_t *parsec_zpotrf_sp = NULL;
    int info = 0;

    parsec_zpotrf_sp = dsparse_zpotrf_sp_New( A, sopalin_data );

    if ( parsec_zpotrf_sp != NULL )
    {
        parsec_enqueue( parsec, (parsec_handle_t*)parsec_zpotrf_sp);
        parsec_context_wait( parsec );
        dsparse_zpotrf_sp_Destruct( parsec_zpotrf_sp );
    }
    return info;
}
