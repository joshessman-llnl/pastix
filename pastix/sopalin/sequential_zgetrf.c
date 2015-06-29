/**
 *
 * @file sequential_zgetrf.c
 *
 *  PaStiX factorization routines
 *  PaStiX is a software package provided by Inria Bordeaux - Sud-Ouest,
 *  LaBRI, University of Bordeaux 1 and IPB.
 *
 * @version 5.1.0
 * @author Pascal Henon
 * @author Xavier Lacoste
 * @author Pierre Ramet
 * @author Mathieu Faverge
 * @date 2013-06-24
 *
 * @precisions normal z -> s d c
 *
 **/
#include "common.h"
#include "sopalin_data.h"
#include "pastix_zcores.h"

void
pastix_static_zgetrf( sopalin_data_t *sopalin_data )
{
    SolverMatrix *datacode = sopalin_data->solvmtx;
    SolverCblk *cblk;
    pastix_int_t  i, ii;
    pastix_int_t tasknbr, *tasktab;
    Task *t;
    FILE *stream = fopen("facto.txt", "w");

    tasknbr = datacode->ttsknbr[0];
    tasktab = datacode->ttsktab[0];

    for (ii=0; ii<tasknbr; ii++){
        i = tasktab[ii];
        t = datacode->tasktab + i;
        cblk = datacode->cblktab + t->cblknum;

        /* Compute task */
        switch( t->taskid )
        {
        case COMP_1D:
            /* Compute */
            core_zgetrfsp1d( datacode, cblk, sopalin_data->diagthreshold );
            break;

        default:
            errorPrint("Taskid unknown for task %ld\n", (long)i);
            EXIT(MOD_SOPALIN,INTERNAL_ERR);
        }
    }
    fclose(stream);
}