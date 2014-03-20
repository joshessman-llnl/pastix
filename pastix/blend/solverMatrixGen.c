#include <stdio.h>
#include <string.h>
#include <strings.h>
#include <math.h>
#include <assert.h>
#include <sys/stat.h>

#include "common.h"
#include "dof.h"
#include "cost.h"
#include "ftgt.h"
#include "symbol.h"
#include "queue.h"
#include "bulles.h"
#include "csc.h"
#include "updown.h"
#include "solver.h"
#include "extendVector.h"
#include "elimin.h"
#include "cand.h"
#include "blendctrl.h"
#include "simu.h"
#include "solver_check.h"
#include "task.h"
#include "fanboth2.h"
#include "solverRealloc.h"
#include "solver_io.h"
#include "solverMatrixGen.h"

/*#define DEBUG_PRIO*/

void working_array_boundaries(SolverMatrix *solvmtx, MPI_Comm pastix_comm);

void printSymbolMatrix(FILE *file, SymbolMatrix *symbptr)
{
    pastix_int_t i, j;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        fprintf(file, "CBLK %ld [%ld : %ld ] \n",(long)i, (long)symbptr->cblktab[i].fcolnum, (long)symbptr->cblktab[i].lcolnum);
        for(j=symbptr->cblktab[i].bloknum;j<symbptr->cblktab[i+1].bloknum;j++)
            fprintf(file, "--BLOK %ld [%ld : %ld ]\n", (long)j, (long)symbptr->bloktab[j].frownum, (long)symbptr->bloktab[j].lrownum);
        fprintf(file, "\n");
    }
}

void build_smx(UpDownVector          *updovct,
               const SymbolMatrix    *symbptr,
               const SimuCtrl        *simuptr,
               const BlendCtrl *const ctrl,
               const Dof       *const dofptr)
{
    pastix_int_t i, j;
    pastix_int_t cursor = 0;
    pastix_int_t xcolnbr = 0;
    pastix_int_t Xnbr = 0;
    pastix_int_t localXnbr = 0;
    pastix_int_t delta;

    for(i=0;i<symbptr->cblknbr;i++)
    {
#ifdef DOF_CONSTANT
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;
        /*if(cbprtab[i] == ctrl->procnum)*/
        if(simuptr->bloktab[symbptr->cblktab[i].bloknum].ownerclust == ctrl->clustnum)
            localXnbr += delta;
        Xnbr += delta;
#else
        EXIT(MOD_BLEND,INTERNAL_ERR);
#endif
    }

    /** We build the whole second member **/
    for(i=0;i<symbptr->cblknbr;i++)
    {
        /** Compute xcolnbr **/
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;


        /** Add delta in the ODB variable **/
        xcolnbr = 0;
        for(j=symbptr->cblktab[i].bloknum+1;j<symbptr->cblktab[i+1].bloknum;j++)
        {
            xcolnbr += (symbptr->bloktab[j].lrownum - symbptr->bloktab[j].frownum + 1)*dofptr->noddval;
        }
        /** We only count non diagonal terms **/
        /** Now add the height of the cblk (-1 for diagonal term) in the DIAGONAL variables **/
        xcolnbr += delta-1;
    }

    /** Now we fill the local second member **/
    updovct->sm2xsze = localXnbr;
    updovct->sm2xnbr = 1;
    updovct->sm2xtab = NULL;

    /* Find the sm2xmax = cblk the broadest on other proc */
    updovct->sm2xmax = 0;
    for(i=0;i<symbptr->cblknbr;i++)
    {
        delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;
        if(updovct->sm2xmax < delta)
            updovct->sm2xmax = delta;
    }

    j = 0;
    for(i=0;i<symbptr->cblknbr;i++)
        /*if(cbprtab[i] == ctrl->procnum)*/
        if(simuptr->bloktab[symbptr->cblktab[i].bloknum].ownerclust == ctrl->clustnum)
        {
            delta = (symbptr->cblktab[i].lcolnum - symbptr->cblktab[i].fcolnum + 1)*dofptr->noddval;

            updovct->cblktab[j].sm2xind = cursor;
            j++;
            cursor += delta;
        }
}



pastix_int_t *
solverMatrixGen(const pastix_int_t clustnum,
                              SolverMatrix *solvmtx,
                              const SymbolMatrix *symbmtx,
                              const SimuCtrl * simuctrl,
                              const BlendCtrl * ctrl,
                              const Dof * dofptr)
{
    pastix_int_t            p, c;
    pastix_int_t            cursor, cursor2;
    pastix_int_t            flag;
    pastix_int_t            i, j, jloc, k;
    pastix_int_t            ftgtnum          = 0;
    pastix_int_t            coefnbr          = 0;
    pastix_int_t            coefind          = 0;
    pastix_int_t            nodenbr          = 0;
    pastix_int_t            odb_nbr          = 0;
    pastix_int_t            cblknum          = 0;
    pastix_int_t            bloknum          = 0;
    pastix_int_t            tasknum          = 0;
    pastix_int_t            facebloknum      = 0;
    pastix_int_t            delta            = 0;
    pastix_int_t            indnbr           = 0;
    pastix_int_t          * cblklocalnum     = NULL;
    pastix_int_t          * bloklocalnum     = NULL;
    pastix_int_t          * tasklocalnum     = NULL;
    pastix_int_t          * ftgtlocalnum     = NULL;
    pastix_int_t          * bcofind          = NULL;
    pastix_int_t          * clust_mask       = NULL;
    pastix_int_t          * clust_first_cblk = NULL;
    pastix_int_t          * clust_highest    = NULL;
    pastix_int_t          * uprecvcblk       = NULL;
    pastix_int_t            flaglocal        = 0;

    /** Set procnum and procnbr **/
    /*solvmtx->procnum = procnum;
     solvmtx->procnbr = ctrl->procnbr;*/
#ifdef PASTIX_DYNSCHED
    solvmtx->btree    = ctrl->btree;
#endif
    solvmtx->clustnum = ctrl->clustnum;
    solvmtx->clustnbr = ctrl->clustnbr;
    solvmtx->procnbr  = ctrl->total_nbcores;
    solvmtx->thrdnbr  = ctrl->local_nbthrds;
    solvmtx->bublnbr  = ctrl->local_nbctxts;
    solvmtx->ftgtcnt  = simuctrl->ftgtcnt;
#ifdef STARPU_GET_TASK_CTX
    solvmtx->starpu_subtree_nbr = symbmtx->starpu_subtree_nbr;
#endif

    if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
    {
        fprintf(stdout, "NUMBER of THREAD %ld \n", (long) solvmtx->thrdnbr );
        fprintf(stdout, "NUMBER of BUBBLE %ld \n", (long) solvmtx->bublnbr );
    }

    /* Copy the vector used to get a cluster number from a processor number */
    MALLOC_INTERN(solvmtx->proc2clust, solvmtx->procnbr, pastix_int_t);
    memcpy(solvmtx->proc2clust, ctrl->core2clust, sizeof(pastix_int_t)*solvmtx->procnbr);

    /** Be sure initialized **/
    solvmtx->cpftmax = 0;
    solvmtx->bpftmax = 0;
    solvmtx->coefmax = 0;

    /*
     * Compute local indices to compress the symbol information into solver
     */
    {
        pastix_int_t *localindex;

        MALLOC_INTERN(localindex, ctrl->clustnbr, pastix_int_t);
        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );

        /* Compute local number of tasks on each cluster */
        MALLOC_INTERN(tasklocalnum, simuctrl->tasknbr, pastix_int_t);
        for(i=0; i<simuctrl->tasknbr; i++) {
            c = simuctrl->bloktab[simuctrl->tasktab[i].bloknum].ownerclust;

            tasklocalnum[i] = localindex[c];
            localindex[c]++;
        }
        solvmtx->tasknbr = localindex[clustnum];

        /* Compute the local numbering of the bloks and cblks on each cluster */
        MALLOC_INTERN(bloklocalnum, symbmtx->bloknbr, pastix_int_t);
        MALLOC_INTERN(cblklocalnum, symbmtx->cblknbr, pastix_int_t);

        memset( localindex, 0, ctrl->clustnbr * sizeof(pastix_int_t) );
        cblknum = 0;
        for(i=0; i<symbmtx->cblknbr; i++)
        {
            flaglocal = 0;
            for(j = symbmtx->cblktab[i].bloknum; j<symbmtx->cblktab[i+1].bloknum; j++)
            {
                c = simuctrl->bloktab[j].ownerclust;
                bloklocalnum[j] = localindex[c];
                localindex[c]++;

                if (c == clustnum)
                    flaglocal = 1;
            }

            if(flaglocal) {
                cblklocalnum[i] = cblknum;
                cblknum++;
            }
            else {
                cblklocalnum[i] = -1;
            }
        }
        solvmtx->bloknbr = localindex[clustnum];
        solvmtx->cblknbr = cblknum;

        memFree_null(localindex);
    }

    /*
     * Fill in bloktab and cblktab
     */

    /* Allocate the cblktab and bloktab with the computed size */
    MALLOC_INTERN(solvmtx->cblktab, solvmtx->cblknbr+1, SolverCblk);
    MALLOC_INTERN(solvmtx->bloktab, solvmtx->bloknbr,   SolverBlok);
    {
        SolverCblk *solvcblk = solvmtx->cblktab;
        SolverBlok *solvblok = solvmtx->bloktab;
        SymbolCblk *symbcblk = symbmtx->cblktab;
        SymbolBlok *symbblok = symbmtx->bloktab;
        SimuBlok   *simublok = simuctrl->bloktab;

        cblknum = 0;
        bloknum = 0;
        nodenbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++, symbcblk++)
        {
            pastix_int_t fbloknum = symbcblk[0].bloknum;
            pastix_int_t lbloknum = symbcblk[1].bloknum;

            flaglocal = 0;
            cursor = bloknum;

            for( j=fbloknum; j<lbloknum; j++, symbblok++, simublok++ ) {

                if(simublok->ownerclust == clustnum)
                {
                    flaglocal = 1;
                    solvblok->frownum = symbblok->frownum * dofptr->noddval;
                    solvblok->lrownum = symbblok->lrownum * dofptr->noddval + dofptr->noddval - 1;
                    solvblok->cblknum = cblklocalnum[symbblok->cblknum];
                    bloknum ++; solvblok++;
                }
            }
            if(flaglocal)
            {
                solvcblk->fcolnum  = symbcblk->fcolnum * dofptr->noddval;
                solvcblk->lcolnum  = symbcblk->lcolnum * dofptr->noddval + dofptr->noddval - 1;
                solvcblk->bloknum  = cursor;
                solvcblk->coeftab  = NULL;
                solvcblk->ucoeftab = NULL;
                nodenbr += symbcblk->lcolnum - symbcblk->fcolnum + 1;
                cblknum++; solvcblk++;
            }
        }

        solvmtx->nodenbr = nodenbr;

        assert( solvmtx->cblknbr == cblknum );
        assert( solvmtx->bloknbr == bloknum );

        if (cblknum > 0)
        {
            /*  virtual cblk to avoid side effect in the loops on cblk bloks */
            solvcblk->fcolnum  = solvcblk->lcolnum + 1;
            solvcblk->lcolnum  = solvcblk->lcolnum + 1;
            solvcblk->bloknum  = bloknum;
            solvcblk->coeftab  = NULL;
            solvcblk->ucoeftab = NULL;
        }
    }

    /*
     * Fill in tasktab
     */
    MALLOC_INTERN(solvmtx->tasktab, solvmtx->tasknbr+1, Task);
    {
        SimuTask *simutask = simuctrl->tasktab;
        Task     *solvtask = solvmtx->tasktab;

        tasknum = 0;
        ftgtnum = 0;
        indnbr  = 0;

        for(i=0; i<simuctrl->tasknbr; i++, simutask++)
        {
            if( simuctrl->bloktab[ simutask->bloknum ].ownerclust == clustnum )
            {
                assert( tasknum == tasklocalnum[i] );

                solvtask->taskid  = simutask->taskid;
                solvtask->prionum = simutask->prionum;
                solvtask->cblknum = cblklocalnum[ simutask->cblknum ];
                solvtask->bloknum = bloklocalnum[ simutask->bloknum ];
                solvtask->ftgtcnt = simutask->ftgtcnt;
                solvtask->ctrbcnt = simutask->ctrbcnt;
                solvtask->indnum  = indnbr;

                /*
                 * Count number of index needed in indtab:
                 *  => number of off-diagonal block below the block (included the block itself)
                 */
                odb_nbr = symbmtx->cblktab[ simutask->cblknum + 1 ].bloknum - simutask->bloknum - 1;

                switch(solvtask->taskid)
                {
                case COMP_1D:
                    indnbr += (odb_nbr*(odb_nbr+1))/2;
                    break;
                default:
                    fprintf(stderr, "solverMatrixgen: Error no task type \n");
                    EXIT(MOD_BLEND,INTERNAL_ERR);
                }

                tasknum++; solvtask++;
            }
        }
        assert(tasknum == solvmtx->tasknbr);

        /* One more task to avoid side effect */
        solvtask->taskid  = -1;
        solvtask->prionum = -1;
        solvtask->cblknum = solvmtx->cblknbr+1;
        solvtask->bloknum = solvmtx->bloknbr+1;
        solvtask->ftgtcnt = 0;
        solvtask->ctrbcnt = 0;
        solvtask->indnum  = indnbr;
    }

    /********************************************/
    /* Fill the processor tasktab indice vector */
    /********************************************/
    /* Number of processor in this cluster */
    k = solvmtx->bublnbr;
    MALLOC_INTERN(solvmtx->ttsknbr, k, pastix_int_t);

    for(p = 0;p<k;p++)
    {
        solvmtx->ttsknbr[p] = extendint_Size(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab);
        print_debug(DBG_BUBBLES, "La bulle %d execute %d taches\n", (int)p, (int)solvmtx->ttsknbr[p]);
    }

    MALLOC_INTERN(solvmtx->ttsktab, k, pastix_int_t *);
    for(p = 0;p<k;p++)
    {
#ifdef PASTIX_DYNSCHED
        pastix_int_t min = INTVALMAX;
        pastix_int_t max = 0;
#endif

        if(solvmtx->ttsknbr[p] > 0)
        {
            MALLOC_INTERN(solvmtx->ttsktab[p], solvmtx->ttsknbr[p], pastix_int_t);
        }
        else
            solvmtx->ttsktab[p] = NULL;

        for(i=0;i<solvmtx->ttsknbr[p];i++)
        {
            j    = extendint_Read(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab, i);
            jloc = tasklocalnum[j];
            solvmtx->ttsktab[p][i] = jloc;

#if (defined PASTIX_DYNSCHED) || (defined TRACE_SOPALIN)
            solvmtx->tasktab[jloc].threadid = p;
#endif
#ifdef TRACE_SOPALIN
            solvmtx->tasktab[jloc].fcandnum = ctrl->candtab[simuctrl->tasktab[j].cblknum].fcandnum;
            solvmtx->tasktab[jloc].lcandnum = ctrl->candtab[simuctrl->tasktab[j].cblknum].lcandnum;
            solvmtx->tasktab[jloc].id       = simuctrl->tasktab[j].cblknum;
#endif

#ifdef PASTIX_DYNSCHED
            if ( solvmtx->tasktab[jloc].prionum > max )
                max = solvmtx->tasktab[jloc].prionum;
            if ( solvmtx->tasktab[jloc].prionum < min )
                min = solvmtx->tasktab[jloc].prionum;
#endif
        }

#ifdef PASTIX_DYNSCHED
        solvmtx->btree->nodetab[p].priomin = min;
        solvmtx->btree->nodetab[p].priomax = max;
#endif
    }

#ifdef PRINT_ORDOTASK
    {
        FILE *ordofile;
        char  ordofilename[250];
        sprintf(ordofilename, "Ordo.%02d", clustnum);
        ordofile = fopen(ordofilename, "w");
        for(p = 0;p<k;p++)
            for(i=0;i<solvmtx->ttsknbr[p];i++)
            {
                j = extendint_Read(simuctrl->proctab[simuctrl->clustab[clustnum].fprocnum + p].tasktab, i);
                fprintf(ordofile, "%ld %ld\n", j, tasklocalnum[j]);
            }
        fclose(ordofile);
    }
#endif

    /*******************/
    /** Fill ftgttab  **/
    /*******************/

    solvmtx->ftgtnbr = 0;
    for(c=0;c < ctrl->clustnbr;c++)
    {
        if(c != clustnum)
            solvmtx->ftgtnbr += extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));
    }

    if(solvmtx->ftgtnbr > 0)
    {
        MALLOC_INTERN(solvmtx->ftgttab, solvmtx->ftgtnbr, FanInTarget);
    }
    else
        solvmtx->ftgttab = NULL;
    MALLOC_INTERN(ftgtlocalnum, simuctrl->bloktab[symbmtx->bloknbr].ftgtnum, pastix_int_t);
    memset(ftgtlocalnum, -1, sizeof(pastix_int_t)*simuctrl->bloktab[symbmtx->bloknbr].ftgtnum);
    cursor = 0;
    for(c=0;c<ctrl->clustnbr;c++)
        for(i=0;i<extendint_Size(&(simuctrl->clustab[clustnum].ftgtsend[c]));i++)
        {
            ftgtnum = extendint_Read(&(simuctrl->clustab[clustnum].ftgtsend[c]), i);
            ftgtlocalnum[ftgtnum] = cursor;
            memcpy(solvmtx->ftgttab[cursor].infotab, simuctrl->ftgttab[ftgtnum].ftgt.infotab, MAXINFO*sizeof(pastix_int_t));


#ifdef DOF_CONSTANT
            solvmtx->ftgttab[cursor].infotab[FTGT_FCOLNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LCOLNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LCOLNUM] += dofptr->noddval - 1;
            solvmtx->ftgttab[cursor].infotab[FTGT_FROWNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LROWNUM] *= dofptr->noddval;
            solvmtx->ftgttab[cursor].infotab[FTGT_LROWNUM] += dofptr->noddval - 1;
#endif

            solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST] = tasklocalnum[solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST]];

            solvmtx->ftgttab[cursor].infotab[FTGT_BLOKDST] = bloklocalnum[solvmtx->ftgttab[cursor].infotab[FTGT_BLOKDST]];
            /* Reinit ftgt ctrbcnt */
            solvmtx->ftgttab[cursor].infotab[FTGT_CTRBCNT] = solvmtx->ftgttab[cursor].infotab[FTGT_CTRBNBR];
            /** Allocation for FanInTarget not assured by blend (done by sopalin)**/
            solvmtx->ftgttab[cursor].coeftab = NULL;

            /*if(p == 2)
             fprintf(stdout, "Ftgt %ld prio %ld to task %ld \n", (long)cursor, (long)solvmtx->ftgttab[cursor].infotab[FTGT_PRIONUM],
             (long)solvmtx->ftgttab[cursor].infotab[FTGT_TASKDST]);*/


            cursor++;
        }


    /*********************/
    /*    Fill indtab    */
    /*********************/
    solvmtx->indnbr   = indnbr;
    solvmtx->indtab   = NULL;
    bcofind           = NULL;
    if (indnbr)
        MALLOC_INTERN(solvmtx->indtab, indnbr, pastix_int_t);

    tasknum     = 0;
    indnbr      = 0;

    for(i=0;i<simuctrl->tasknbr;i++)
    {
        if(simuctrl->bloktab[simuctrl->tasktab[i].bloknum].ownerclust == clustnum)
        {
            ASSERTDBG(tasklocalnum[i] == tasknum,MOD_BLEND);
            ASSERTDBG(indnbr == solvmtx->tasktab[tasklocalnum[i]].indnum, MOD_BLEND);
            ASSERTDBG(bloklocalnum[simuctrl->tasktab[i].bloknum] == solvmtx->tasktab[tasklocalnum[i]].bloknum,MOD_BLEND);
            ASSERTDBG(cblklocalnum[simuctrl->tasktab[i].cblknum] == solvmtx->tasktab[tasklocalnum[i]].cblknum,MOD_BLEND);

            switch(simuctrl->tasktab[i].taskid)
            {
            case COMP_1D:
                for(bloknum = symbmtx->cblktab[simuctrl->tasktab[i].cblknum].bloknum+1;
                    bloknum < symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum; bloknum++)
                {
                    facebloknum = 0;
                    for(j=bloknum;j<symbmtx->cblktab[simuctrl->tasktab[i].cblknum+1].bloknum;j++)
                    {
                        facebloknum = symbolGetFacingBloknum(symbmtx, bloknum, j, facebloknum, ctrl->ricar);
                        /*#ifdef NAPA*/
                        if(facebloknum >= 0)
                        {
                            /*#endif*/
                            if(simuctrl->bloktab[facebloknum].ownerclust!=clustnum)
                            {
                                solvmtx->indtab[indnbr] = ftgtlocalnum[CLUST2INDEX(facebloknum, clustnum)];
#ifdef DEBUG_PRIO
                                solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]
                                    = solvmtx->tasktab[tasklocalnum[i]].prionum;
                                /*fprintf(stdout, "SOLV Task1D %ld FTGT %ld  taskprio %ld ftgtprio %ld \n", (long)tasklocalnum[i], (long)solvmtx->indtab[indnbr] , (long)solvmtx->tasktab[tasklocalnum[i]].prionum, (long)solvmtx->ftgttab[solvmtx->indtab[indnbr]].infotab[FTGT_PRIONUM]);*/
#endif

                                ASSERTDBG(solvmtx->indtab[indnbr] < solvmtx->ftgtnbr,MOD_BLEND);
                            }
                            else
                            {
                                solvmtx->indtab[indnbr] = -tasklocalnum[simuctrl->bloktab[facebloknum].tasknum];

#ifdef DEBUG_BLEND
                                if(!(-solvmtx->indtab[indnbr] < solvmtx->tasknbr))
                                    fprintf(stderr, "facetasknum %ld tasklocalnum %ld \n",
                                            (long)simuctrl->bloktab[facebloknum].tasknum, (long)-solvmtx->indtab[indnbr]);
                                ASSERT(-solvmtx->indtab[indnbr] < solvmtx->tasknbr,MOD_BLEND);
#endif
                            }
                            indnbr++;
                            ASSERTDBG(indnbr <= solvmtx->indnbr,MOD_BLEND);
                            /*#ifdef NAPA*/
                        }
                        else
                        { /** THE FACING BLOCK DO NOT EXIST **/
                            solvmtx->indtab[indnbr] =  solvmtx->ftgtnbr+1;
                            indnbr++;
                        }
                        /*#endif*/
                    }
                }
                break;
            default:
                fprintf(stderr, "Error in solverMatrixgen \n");
                EXIT(MOD_BLEND,INTERNAL_ERR);
            }
            tasknum++;
        }
    }
    ASSERTDBG(indnbr  == solvmtx->indnbr,  MOD_BLEND);

    /********************/
    /** Fill Solver    **/
    /** cblk and blok  **/
    /********************/
    cblknum = 0;
    bloknum = 0;
    coefnbr = 0;
    coefind = 0;
    for(i=0;i<solvmtx->cblknbr;i++)
    {
        coefind = 0;
        solvmtx->cblktab[i].stride   = 0;
        solvmtx->cblktab[i].procdiag = -1;
        for(j=solvmtx->cblktab[i].bloknum;j<solvmtx->cblktab[i+1].bloknum;j++)
        {
            /* Solvmtx is already expanded in dll */
            delta =  solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum +1;
            solvmtx->cblktab[i].stride += delta;
            solvmtx->bloktab[j].coefind = coefind;
            coefind += delta;
        }
        coefnbr += (solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum +1) * solvmtx->cblktab[i].stride;
        coefind = coefnbr;
    }
    solvmtx->coefnbr = coefnbr;

    /*****************************************/
    /**  Find coefmax, cpftmax and bpftmax  **/
    /*****************************************/
    /** Find bpftmax **/
    /* bpftmax is the number of coef of the largest block target in reception */
    solvmtx->bpftmax = 0;
    /* the largest block target is the largest local block */


    /** Find cpftmax **/
    /* cpftmax is the number of coef of the largest fan in target in reception */
    solvmtx->cpftmax = 0;

    {
        /***** Find coefmax *****
         * coefmax is the number of coef of the largest temporary block used
         * to compute contribution block on the CPUs.
         * It can be seen as the maximum surface of the
         * C matrix in the GEMM operations.
         */
        pastix_int_t max_m = 0;
        pastix_int_t max_n = 0;

        solvmtx->coefmax = 0;

        if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES) {
            for(i=0;i<solvmtx->tasknbr;i++) {
                if(solvmtx->tasktab[i].taskid == COMP_1D) {
                    delta = 0;
                    assert(solvmtx->tasktab[i].cblknum >= 0);
                    for(j=solvmtx->tasktab[i].bloknum;
                        j<solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].bloknum;
                        j++) {
                        coefind = solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum+1;
#ifdef PASTIX_ESC
                        while(((j+1) < solvmtx->cblktab[solvmtx->tasktab[i].cblknum+1].bloknum)
                              && (solvmtx->bloktab[j].cblknum == solvmtx->bloktab[j+1].cblknum)) {
                            j++;
                            coefind += solvmtx->bloktab[j].lrownum - solvmtx->bloktab[j].frownum+1;
                        }
#endif
                        if(coefind > delta)
                            delta = coefind;
                    }
                    coefnbr = solvmtx->cblktab[solvmtx->tasktab[i].cblknum].stride * delta;
                    if(coefnbr > solvmtx->coefmax) {
                        solvmtx->coefmax = coefnbr;
                        max_m = solvmtx->cblktab[solvmtx->tasktab[i].cblknum].stride;
                        max_n = delta;
                    }
                }
            }
            /* Maximum size of diagonal blocks */
            for(i=0;i<solvmtx->cblknbr;i++) {
                delta = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                if(delta*delta > solvmtx->coefmax) {
                    solvmtx->coefmax = delta*delta;
                    max_m = delta;
                    max_n = delta;
                }
            }

            fprintf(stderr, "Actual coefmax = %ld (%ld x %ld)\n",
                    (long)solvmtx->coefmax, (long)max_m, (long)max_n );
        }

        /* First compute the maximum size of contribution block */
        solvmtx->coefmax = 0;
        {
            pastix_int_t itercblk;
            pastix_int_t m, n;
            for (itercblk = 0; itercblk < solvmtx->cblknbr; itercblk++) {
                pastix_int_t stride = solvmtx->cblktab[ itercblk ].stride;
                for(i=solvmtx->cblktab[ itercblk ].bloknum+1;
                    i<solvmtx->cblktab[ itercblk + 1].bloknum;i++) {
                    m = stride - solvmtx->bloktab[i].coefind;
                    n = solvmtx->bloktab[i].lrownum - solvmtx->bloktab[i].frownum+1;
#ifdef PASTIX_ESC
                    while(((i+1) < solvmtx->cblktab[itercblk+1].bloknum)
                          && (solvmtx->bloktab[i].cblknum == solvmtx->bloktab[i+1].cblknum)) {
                        i++;
                        n += solvmtx->bloktab[n].lrownum - solvmtx->bloktab[n].frownum+1;
                    }
#endif
                    delta = m * n;
                    if(delta > solvmtx->coefmax) {
                        solvmtx->coefmax = delta;
                        max_m = m;
                        max_n = n;
                    }
                }
                /* kernel_trsm require COLNBR * (stride - COLNBR + 1) in LDLt */
                /* horizontal dimension */
                n = solvmtx->cblktab[itercblk].lcolnum -
                    solvmtx->cblktab[itercblk].fcolnum + 1;
                /* vertical dimension */
                m = stride - n + 1;
                delta = m * n;
                if(delta > solvmtx->coefmax) {
                    solvmtx->coefmax = delta;
                    max_m = m;
                    max_n = n;
                }
            }
        }

        if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES) {
            fprintf(stderr, "New suggested coefmax = %ld (%ld x %ld)\n",
                    (long)solvmtx->coefmax, (long)max_m, (long)max_n );
            /* First compute the maximum size of contribution block */
            {
                pastix_int_t max = 0;
                pastix_int_t n;
                for(i=0;i<solvmtx->cblknbr-1;i++) {
                    n = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                    delta = n * 64;
                    if(delta > max) {
                        max = delta;
                        max_m = n;
                        max_n = 64;
                    }
                }
                fprintf(stderr, "Max diagblock coefmax without shur = %ld (%ld x %ld)\n",
                        (long)max, (long)max_m, (long)max_n );

                n = solvmtx->cblktab[i].lcolnum - solvmtx->cblktab[i].fcolnum+1;
                delta = n * 64;
                if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_YES)
                    fprintf(stderr, "Max diagblock on shur = %ld (%ld x %ld)\n",
                            (long)delta, (long)n, (long)64 );
            }
        }
    }

    /** Find the cpftmax **/
    /* OIMBE on peut trouver le bon : flemmard */
    solvmtx->cpftmax = 0;
    for(i=0;i<simuctrl->ftgtnbr;i++)
    {
        if((simuctrl->ftgttab[i].ftgt.infotab[FTGT_CTRBNBR]>0))
            /*&& (proc2clust[simuctrl->ftgttab[i].ftgt.infotab[FTGT_PROCDST]] == clustnum))*/
        {
            coefnbr = (simuctrl->ftgttab[i].ftgt.infotab[FTGT_LCOLNUM] - simuctrl->ftgttab[i].ftgt.infotab[FTGT_FCOLNUM] + 1)*dofptr->noddval * (simuctrl->ftgttab[i].ftgt.infotab[FTGT_LROWNUM] - simuctrl->ftgttab[i].ftgt.infotab[FTGT_FROWNUM] + 1)*dofptr->noddval;
            if(coefnbr > solvmtx->cpftmax)
                solvmtx->cpftmax = coefnbr;
        }
    }

    /** Find the bpftmax **/
    solvmtx->bpftmax = 0;
    for(i=0;i<simuctrl->tasknbr;i++)
    {
        if(simuctrl->tasktab[i].taskid == E1 || simuctrl->tasktab[i].taskid == E2)
            if(simuctrl->bloktab[simuctrl->tasktab[i].bloknum].ownerclust == clustnum)
            {
                delta = (symbmtx->cblktab[simuctrl->tasktab[i].cblknum].lcolnum - symbmtx->cblktab[simuctrl->tasktab[i].cblknum].fcolnum+1)
                    * dofptr->noddval;
                coefnbr = delta * (symbmtx->bloktab[simuctrl->tasktab[i].bloknum2].lrownum - symbmtx->bloktab[simuctrl->tasktab[i].bloknum2].frownum+1)
                    *dofptr->noddval;
                if(coefnbr > solvmtx->bpftmax)
                    solvmtx->bpftmax = coefnbr;

            }
    }



    /** Find the nbftmax **/
    solvmtx->nbftmax = 0;
    for(i=0;i<simuctrl->tasknbr;i++)
        if(simuctrl->tasktab[i].ftgtcnt> solvmtx->nbftmax)
            solvmtx->nbftmax = simuctrl->tasktab[i].ftgtcnt;

    /** Find the area max **/

    /** We search the biggest sum of AUB that can be send to the local cluster from a cblk on another cluster **/
    /* @@@@@@@@@@@@ TODO  ****/

    /** We search the biggest cblk that can receive some ftgt **/
    solvmtx->arftmax = 0;
    for(i=0;i<simuctrl->cblknbr;i++)
    {
        pastix_int_t size = 0;
        delta = (symbmtx->cblktab[i].lcolnum -  symbmtx->cblktab[i].fcolnum + 1)*dofptr->noddval;
        if(ctrl->candtab[i].fccandnum != ctrl->candtab[i].lccandnum)
        {
            for(j=symbmtx->cblktab[i].bloknum;j<symbmtx->cblktab[i+1].bloknum;j++)
                size += delta * (symbmtx->bloktab[j].lrownum-symbmtx->bloktab[j].frownum+1) * dofptr->noddval;

            if(size> solvmtx->arftmax)
                solvmtx->arftmax = size;
        }
    }


    if (ctrl->iparm[IPARM_VERBOSE]>API_VERBOSE_NO)
        fprintf(stdout, "COEFMAX %ld CPFTMAX %ld BPFTMAX %ld NBFTMAX %ld ARFTMAX %ld \n", (long)solvmtx->coefmax, (long)solvmtx->cpftmax,
                (long)solvmtx->bpftmax, (long)solvmtx->nbftmax, (long)solvmtx->arftmax);


    /****************************************/
    /** Compute the information for the    **/
    /** Forward and Backward triangular    **/
    /** Solution                           **/
    /****************************************/

    /* Pour l'instant uniquement si on est en 1d */
    if (ctrl->level2D == 0)
    {

        /** The initial symbol matrix is not expanded **/
        nodenbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            nodenbr += (symbmtx->cblktab[i].lcolnum-symbmtx->cblktab[i].fcolnum+1)*dofptr->noddval;
        solvmtx->updovct.gnodenbr= nodenbr;
        /*fprintf(stderr," GNODENBR %ld \n", (long)solvmtx->updovct.gnodenbr);*/

        /** Build the browtabs for each diagonal block **/
        MALLOC_INTERN(solvmtx->updovct.cblktab, solvmtx->cblknbr,UpDownCblk);
        cursor = 0;
        MALLOC_INTERN(clust_mask,       ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_first_cblk, ctrl->clustnbr, pastix_int_t);
        MALLOC_INTERN(clust_highest,    ctrl->clustnbr, pastix_int_t);

        solvmtx->updovct.downmsgnbr = 0;

        for(i=0;i<symbmtx->cblknbr;i++)
        {
            pastix_int_t brownbr;

            /*if(cbprtab[i] == clustnum)*/
            if(simuctrl->bloktab[symbmtx->cblktab[i].bloknum].ownerclust == clustnum)
            {
                /*** Compute the list of clusters in the BROW (each cluster is list one time) ***/
                bzero(clust_mask, sizeof(pastix_int_t)*ctrl->clustnbr);
                bzero(clust_first_cblk, sizeof(pastix_int_t)*ctrl->clustnbr);
                for(j=0;j<ctrl->clustnbr;j++)
                    /*clust_highest[j] = - simuctrl->cblknbr;*/
                    clust_highest[j] = -1;


                brownbr = 0;
                for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                {
                    pastix_int_t cluster;
                    pastix_int_t cblk;
                    cluster = simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust;
                    cblk = ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]];
#ifdef DEBUG_M
                    ASSERT( ctrl->candtab[cblk].treelevel   <= 0,   MOD_BLEND);
                    ASSERT( ctrl->candtab[cblk].distrib     == D1,  MOD_BLEND);
                    ASSERT( simuctrl->tasktab[cblk].cblknum == cblk,MOD_BLEND);
#endif

                    /*if( ctrl->candtab[cblk].treelevel >= clust_highest[cluster] )*/
                    if( simuctrl->tasktab[cblk].prionum >= clust_highest[cluster] )
                    {
                        /*clust_highest[cluster] = ctrl->candtab[cblk].treelevel;*/
                        clust_highest[cluster] = simuctrl->tasktab[cblk].prionum;
                        clust_first_cblk[cluster] = cblk;
                    }

                    if(clust_mask[cluster] == 0)
                    {
                        clust_mask[cluster] = 1;
                        brownbr++;
                    }
                }

                solvmtx->updovct.cblktab[cursor].browprocnbr = brownbr;
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browproctab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browproctab = NULL;

                if(clust_mask[clustnum] == 1)
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr-1;
                else
                    solvmtx->updovct.cblktab[cursor].msgnbr = brownbr;
                solvmtx->updovct.downmsgnbr += solvmtx->updovct.cblktab[cursor].msgnbr;

                /*** Alloc the vector that will contain the global cblknum with the max priority for each processor in browproctab  ***/
                if(solvmtx->updovct.cblktab[cursor].browprocnbr>0)
                {
                    MALLOC_INTERN(solvmtx->updovct.cblktab[cursor].browcblktab,
                                  solvmtx->updovct.cblktab[cursor].browprocnbr,
                                  pastix_int_t);
                }
                else
                    solvmtx->updovct.cblktab[cursor].browcblktab = NULL;

                brownbr = 0;
                for(j=0;j<ctrl->clustnbr;j++)
                    if(clust_mask[j] == 1)
                    {
                        solvmtx->updovct.cblktab[cursor].browproctab[brownbr]   = j;
                        solvmtx->updovct.cblktab[cursor].browcblktab[brownbr++] = clust_first_cblk[j];
                    }

                solvmtx->updovct.cblktab[cursor].ctrbnbr = (pastix_int_t)ctrl->egraph->verttab[i].innbr;
                cursor++;
            }
        }

        /********************************************************************/
        /*** Find the list of local blocks in front of a diagonal blocks ****/
        /********************************************************************/
        cursor  = 0;
        cursor2 = 0;
        MALLOC_INTERN(solvmtx->updovct.gcblk2list, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            solvmtx->updovct.gcblk2list[i] = -1;
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                if( simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        flag = 1;
                        solvmtx->updovct.gcblk2list[i] = cursor;
                        cursor++;
                    }
                    cursor2++;
                }
        }
        solvmtx->updovct.gcblk2listnbr = symbmtx->cblknbr;

        MALLOC_INTERN(solvmtx->updovct.listptr, cursor+1, pastix_int_t);
        solvmtx->updovct.listptrnbr = cursor+1;
        MALLOC_INTERN(solvmtx->updovct.listblok, cursor2, pastix_int_t);
        MALLOC_INTERN(solvmtx->updovct.listcblk, cursor2, pastix_int_t);

        cursor  = 0;
        cursor2 = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
        {
            flag = 0;
            for(j=0; j<ctrl->egraph->verttab[i].innbr;j++)
                if( simuctrl->bloktab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]].ownerclust == clustnum)
                {
                    if(flag == 0)
                    {
                        solvmtx->updovct.listptr[cursor2] = cursor;
                        cursor2++;
                        flag = 1;
                    }
#ifdef OOC
                    {
                        pastix_int_t tmp1,tmp2, tmp3;
                        pastix_int_t iter;
                        tmp1 = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                        tmp2 = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
                        for (iter = solvmtx->updovct.listptr[cursor2-1]; iter < cursor; iter ++)
                            if (solvmtx->tasktab[tmp2].prionum <
                                solvmtx->tasktab[solvmtx->updovct.listcblk[iter]].prionum )
                            {
                                /* No problem with using solvmtx->updovct.listcblk[iter]
                                 * during first loop, cursor = 0, we don't use it,
                                 * and we set first.
                                 */
                                tmp3 = solvmtx->updovct.listcblk[iter];
                                solvmtx->updovct.listcblk[iter] = tmp2;
                                tmp2 = tmp3;

                                tmp3 = solvmtx->updovct.listblok[iter];
                                solvmtx->updovct.listblok[iter] = tmp1;
                                tmp1 = tmp3;
                            }
                        solvmtx->updovct.listblok[cursor] = tmp1;
                        solvmtx->updovct.listcblk[cursor] = tmp2;

                    }
#else
                    solvmtx->updovct.listblok[cursor] = bloklocalnum[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j] ];
                    solvmtx->updovct.listcblk[cursor] = cblklocalnum[ctrl->egraph->ownetab[ctrl->egraph->inbltab[ctrl->egraph->verttab[i].innum+j]]];
#endif
                    cursor++;
                }

        }
        solvmtx->updovct.listptr[cursor2] = cursor;
        solvmtx->updovct.listnbr = cursor;

        solvmtx->updovct.loc2globnbr = solvmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.loc2glob, solvmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            if(cblklocalnum[i] >= 0)
                solvmtx->updovct.loc2glob[cblklocalnum[i]] = i;

        memFree_null(clust_mask);
        memFree_null(clust_first_cblk);
        memFree_null(clust_highest);

        /***** Fill lblk2gcblk ******/
        solvmtx->updovct.gcblknbr = symbmtx->cblknbr;
        MALLOC_INTERN(solvmtx->updovct.lblk2gcblk, symbmtx->bloknbr, pastix_int_t);
        for(i=0;i<symbmtx->bloknbr;i++)
            if(simuctrl->bloktab[i].ownerclust == clustnum)
                solvmtx->updovct.lblk2gcblk[bloklocalnum[i]] = symbmtx->bloktab[i].cblknum;

        /* Calcul du nombre de messages a recevoir lors de la remont�e */
        MALLOC_INTERN(uprecvcblk, symbmtx->cblknbr, pastix_int_t);
        for(i=0;i<symbmtx->cblknbr;i++)
            uprecvcblk[i] = 0;
        for (i=0; i<solvmtx->bublnbr; i++)
            for (j=0; j < solvmtx->ttsknbr[i]; j++)
            {
                cblknum = solvmtx->tasktab[solvmtx->ttsktab[i][j]].cblknum;
                for (k =  solvmtx->cblktab[cblknum+1].bloknum-1;
                     k >= solvmtx->cblktab[cblknum].bloknum+1; k--)
                    /* if the contribution is not local */
                    if (solvmtx->bloktab[k].cblknum <= 0)
                        uprecvcblk[solvmtx->updovct.lblk2gcblk[k]] = 1;
            }
        solvmtx->updovct.upmsgnbr = 0;
        for(i=0;i<symbmtx->cblknbr;i++)
            solvmtx->updovct.upmsgnbr += uprecvcblk[i];
        memFree_null(uprecvcblk);

        /*********************************/
        /*     Temporaire                */
        /*  Pour tester descente remonte */
        /*********************************/
        build_smx(&(solvmtx->updovct), symbmtx, simuctrl, ctrl, dofptr);

    }
    /*********************** END TRIANGULAR INFO BUILDING ******************************************/

    memFree_null(cblklocalnum);
    memFree_null(bloklocalnum);
    memFree_null(tasklocalnum);

    if(ftgtlocalnum != NULL)
        memFree_null(ftgtlocalnum);

    return bcofind;
}
