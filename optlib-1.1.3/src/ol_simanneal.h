//  OptLib version 1.0
//    Copyright (C) 2012  David A. Joiner
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef SIMANNEAL
#define SIMANNEAL

#include <stdio.h>
#include <stdlib.h>
#include "ol_mem.h"
#include "ol_rand_tools.h"
#include <math.h>
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

#define SIMANNEAL_VERBOSE 0
#define SIMANNEAL_STEP_METHOD_1 1
#define SIMANNEAL_CHECKPOINT 0
#define SIMANNEAL_CHECKPOINT_FREQ 1

#define SIMANNEAL_ADAPTIVE_FALSE 0
#define SIMANNEAL_ADAPTIVE_TRUE 1

#ifdef USE_OPTMPI
int ebsa_mpi_rank;
int ebsa_mpi_size;
MPI_Status ebsa_mpi_status;
MPI_Comm ebsa_mpi_comm;
#endif

typedef struct {
    double E;
    double *lowest_x;
    double *temp_x;
    double *t_step,*new_x;
    int *index;
    double lowest_E;
    double T,step;
    int iter,restart_counter;
    int n, n_annealers;
    int sync_counter;
    int simanneal_step_type;
} SimAnneal;

void simanneal_set_adaptive(SimAnneal * theSim, int flag);
void simanneal_adjustSteps(SimAnneal * theSim,int n, int n_annealers,double * x, double * step_factor,
		double (*func)(double * x, void * func_data), void * func_data);

void simanneal_init(SimAnneal ** theSim,int n,int n_annealers,double * x,double * step_factor, double T,
        double step,double step_max, double cooling,
        int ITMAX,int restart, double EPS,
        double (*func)(double * x,void *),void * func_data,
        void output(double * x, void * func_data, const char * mesg));

void simanneal_finalize(SimAnneal * theSim);

int simanneal_step(SimAnneal *,int n, double *x,double * step_factor,
        double *E, double T,double step,
        double (*func)(double *x,void *),void *, void output(double * x,
        void * func_data, const char * mesg));

int simanneal_iterate(SimAnneal * theSim,int n, int n_annealers,double * x,double * step_factor,
        double step_max, double cooling,
        int ITMAX,int restart, double EPS,
        double (*func)(double * x,void *),void *,
        void output(double * x, void * func_data, const char *));
void simanneal(SimAnneal ** pSim,int n, int n_annealers,double * x,double * step_factor,
        double T, double step, double step_max, double cooling,
        int ITMAX,int restart, double EPS, int ADAPTIVE,
        double (*func)(double * x,void *),void *,
        void output(double * x, void * func_data, const char *),
        const char * cp_fname,int RECOVER);
#endif
