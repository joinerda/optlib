#include "optlib.h"
#include "ol_fw_mpi.h"
#include "ol_simanneal.h"
#include "ol_genetic_alg.h"
#include "ol_powell.h"
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

OptLibOpts * OPTLIB_CreateOpts() {
    OptLibOpts * theOpts;
    theOpts = (OptLibOpts *)malloc(sizeof(OptLibOpts));
    theOpts->method = OPTLIB_METHOD_GA;
    theOpts->ga_npools = 3;
    theOpts->ga_npop = 100;
    theOpts->ga_nkeep = 10;
    theOpts->ga_itmax = 100;
    theOpts->ga_itmin = 20;
    theOpts->ga_epsilon = 1.0e-4;
    theOpts->ga_nstep = 1;
    theOpts->ga_stepfactor = NULL;
    theOpts->ga_mutationfactor = 0.05;
    theOpts->ga_dominancefactor = 0.1;
    theOpts->powell_itmax = 10000;
    theOpts->powell_epsilon = 1.0e-3;
    theOpts->mcp_ntrials=100;
    theOpts->mcp_scale=10.0;
    theOpts->sa_nannealers=5;
    theOpts->sa_stepfactor = NULL;
    theOpts->sa_T=10.0;
    theOpts->sa_step=0.1;
    theOpts->sa_stepmax=10.0;
    theOpts->sa_coolingfactor=0.90;
    theOpts->sa_itmax=5000000;
    theOpts->sa_restart=2500000;
    theOpts->sa_epsilon=1.0e-3;
    theOpts->output=NULL;
#ifdef USE_OPTMPI
    theOpts->comm=MPI_COMM_WORLD;
#endif
    return theOpts;
}
void OPTLIB_DestroyOpts(OptLibOpts * theOpts) {
    if(theOpts->ga_stepfactor!=NULL) free(theOpts->ga_stepfactor);
    if(theOpts->sa_stepfactor!=NULL) free(theOpts->sa_stepfactor);
    free(theOpts);
}

void OPTLIB_Minimize(
    int n, double * guess,
    double (*func)(double *, void *),
    void * func_data, OptLibOpts * theOpts_in) {

    int destroyOpts=0;
    OptLibOpts * theOpts;

    if(theOpts_in==NULL) {
        destroyOpts=1;
        theOpts = OPTLIB_CreateOpts();
    } else {
        theOpts = theOpts_in;
    }

    if (theOpts->method==OPTLIB_METHOD_GA) {
        fw_struct * func_struct;
        func_struct = (fw_struct*)malloc(sizeof(fw_struct));
        func_struct->func_data = func_data;
        func_struct->func = func;
        func_struct->n = n;
#ifdef USE_OPTMPI
        func_struct->foreman_works = TRUE;
        func_struct->comm = theOpts->comm;
        MPI_Comm_size(theOpts->comm,&(func_struct->size));
        MPI_Comm_rank(theOpts->comm,&(func_struct->rank));
#else
        func_struct->size=1;
        func_struct->rank=0;
#endif

#ifdef USE_OPTMPI
        if(func_struct->rank!=0) {
            fw_worker_main(func_struct);
        } else {
#endif
            if(theOpts->ga_npools>1) {
                genetic_pools(n,theOpts->ga_npools,
                    theOpts->ga_npop, theOpts->ga_nkeep,
                    theOpts->ga_itmax, theOpts->ga_itmin, theOpts->ga_epsilon,
                    theOpts->ga_nstep, theOpts->ga_stepfactor,
                    theOpts->ga_mutationfactor, theOpts->ga_dominancefactor,
                    guess, theOpts->output, func_struct);
            } else {
                genetic_alg(n,theOpts->ga_npop, theOpts->ga_nkeep,
                    theOpts->ga_itmax, theOpts->ga_itmin, theOpts->ga_epsilon,
                    theOpts->ga_nstep, theOpts->ga_stepfactor,
                    theOpts->ga_mutationfactor, theOpts->ga_dominancefactor,
                    guess, theOpts->output, func_struct);
            }
#ifdef USE_OPTMPI
            fw_shutdown_workers(func_struct->size,func_struct->comm);
        } 
#endif
        free(func_struct);
#ifdef USE_OPTMPI
        MPI_Bcast(guess,n,MPI_DOUBLE,0,theOpts->comm);
#endif
    } else if(theOpts->method==OPTLIB_METHOD_SA) {
#ifdef USE_OPTMPI
        ebsa_mpi_comm=theOpts->comm;
#endif
        int i,j;
        SimAnneal * theSim;
        double * anneal_guess;

        anneal_guess = alloc_darray(n*theOpts->sa_nannealers);
        for(j=0;j<theOpts->sa_nannealers;j++) {
            for(i=0;i<n;i++) {
                anneal_guess[j*n+i]=guess[i];
            }
        }

        simanneal(&theSim,n,theOpts->sa_nannealers,anneal_guess,
            theOpts->sa_stepfactor,theOpts->sa_T,theOpts->sa_step,
            theOpts->sa_stepmax,theOpts->sa_coolingfactor,
            theOpts->sa_itmax,theOpts->sa_restart,theOpts->sa_epsilon,
            SIMANNEAL_ADAPTIVE_TRUE,
            func,func_data,theOpts->output,NULL,0);

        for(j=0;j<n;j++) guess[j] = theSim->lowest_x[j];

        simanneal_finalize(theSim);
        free(anneal_guess);

    } else if(theOpts->method==OPTLIB_METHOD_POWELL) {
        powell(n,guess,func,func_data,theOpts->output,
            theOpts->powell_itmax,
            theOpts->powell_epsilon);
    } else if(theOpts->method==OPTLIB_METHOD_MCP) {
#ifdef USE_OPTMPI
        mcp_mpi_comm = theOpts->comm;
#endif
        mc_powell(theOpts->mcp_ntrials,theOpts->mcp_scale,
            NULL,NULL,NULL,
            n,guess,func,func_data,theOpts->output,
            theOpts->powell_itmax,
            theOpts->powell_epsilon);
    } else {
        printf("OPTLIB_Minimize, method not defined\n");
        exit(0);
    }
    if(destroyOpts) {
        OPTLIB_DestroyOpts(theOpts);
    }
}

