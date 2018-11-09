#ifndef OPTLIB_H
#define OPTLIB_H

#include "ol_rand_tools.h"
#include "ol_fw_mpi.h"
#include "ol_genetic_alg.h"
#include "ol_interp.h"
#include "ol_mem.h"
#include "ol_powell.h"
#include "ol_rand_tools.h"
#include "ol_simanneal.h"
#include "ol_transform.h"
#include "ol_config.h"
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

#define OPTLIB_METHOD_GA 0
#define OPTLIB_METHOD_SA 1
#define OPTLIB_METHOD_MCP 2
#define OPTLIB_METHOD_POWELL 3

/*
void mc_powell(int trials, double scale, double * x_attempts,
        double * x_values,
        double * fvalues, int n,double *x,
        double func(double *,void *),void *,
        void output(double * x, void * func_data, const char *),
        int itmax,double eps);
void powell(int n,double *x,double func(double *,void *),void *,
        void output(double * x, void * func_data, const char *),
        int itmax,double eps);
*/


typedef struct {
    int method;
    int ga_npools;
    int ga_npop;
    int ga_nkeep;
    int ga_itmax;
    int ga_itmin;
    double ga_epsilon;
    int ga_nstep;
    double ga_mutationfactor;
    double ga_dominancefactor;
    int powell_itmax;
    double powell_epsilon;
    int mcp_ntrials;
    double mcp_scale;
    double * ga_stepfactor;
    int sa_nannealers;
    double * sa_stepfactor;
    double sa_T;
    double sa_step;
    double sa_stepmax;
    double sa_coolingfactor;
    int sa_itmax;
    int sa_restart;
    double sa_epsilon;
    void (*output)(double * x, void * func_data, const char *);
#ifdef USE_OPTMPI
    MPI_Comm comm;
#endif
} OptLibOpts;

void OPTLIB_Minimize(
    int n, double * guess,
    double (*func)(double *, void *),
    void * func_data, OptLibOpts * theOpts);

OptLibOpts * OPTLIB_CreateOpts();
void OPTLIB_DestroyOpts(OptLibOpts *);

#endif
