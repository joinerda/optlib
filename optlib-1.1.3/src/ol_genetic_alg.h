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

#ifndef GENETIC_ALG
#define GENETIC_ALG


#include <stdio.h>
#include <math.h>
#include "ol_mem.h"
#include "ol_rand_tools.h"
#include "ol_fw_mpi.h"

#define GA_MPI_VERBOSE 0

void quicksort(int n, double * x, int * index);
void genetic_combination(int n, double * child,
    double * mother, double * father,
    int n_step, double * step_factor,double mutation_factor,
    double dominance_factor, double * range);
void genetic_mutation(int n, double *mutation, double *original,
    int n_step, double * step_factor);
void genetic_mix(int n,int n_pop, int n_keep,
    int itmax, int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor,
    double dominance_factor,
    double ** population,
    void output(double * x, void * func_data, const char * mesg),
    fw_struct * func_struct);
void genetic_alg_serial(int n, int n_pop, int n_keep,
    int itmax, int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor,
    double dominance_factor,
    double * guess,
    void output(double * x, void * func_data, const char * mesg),
    double func(double *, void *), void *);
void genetic_alg(int n, int n_pop, int n_keep,
    int itmax, int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor,
    double dominance_factor,
    double * guess,
    void output(double * x, void * func_data, const char * mesg),
    fw_struct * func_struct);
void genetic_pools_serial(int n, int n_pools, int n_pop, int n_keep,
    int itmax, int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor,
    double dominance_factor,
    double * guess, 
    void output(double * x, void * func_data, const char * mesg),
    double func(double *, void *), void *);
void genetic_pools(int n, int n_pools, int n_pop, int n_keep,
    int itmax, int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor,
    double dominance_factor,
    double * guess, 
    void output(double * x, void * func_data, const char * mesg),
    fw_struct * func_struct);

#endif
