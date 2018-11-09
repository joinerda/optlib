
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




#ifndef POWELL
#define POWELL

#ifdef USE_OPTMPI
#include <mpi.h>

int mcp_mpi_size, mcp_mpi_rank;
MPI_Status mcp_mpi_status;
MPI_Comm mcp_mpi_comm;
#endif

#define POWELL_VERBOSE 0

double powell_directional_deriv(
        int n, int i, double ** u,double * x,double func(double * x,void *),void *);
double powell_directional_deriv2(
        int n, int i, double ** u,double * x,double func(double * x,void *),void *);
double powell_findmin(int n,int i,double ** u,
        double ** P,double func(double * x, void *),void *,
        int itmax, double eps);
void mc_powell(int trials, double scale, double * x_attempts,
        double * x_values,
        double * fvalues, int n,double *x,
        double func(double *,void *),void *,
        void output(double * x, void * func_data, const char *),
        int itmax,double eps);
void powell(int n,double *x,double func(double *,void *),void *,
        void output(double * x, void * func_data, const char *),
        int itmax,double eps);
double vector_distance(int n, double * x, double * y);
double vector_length(int n, double * x);
double powell_newton_min(int n, int i, double ** u, double ** P,
        double func(double *x,void *),void *,int itmax, double eps,int * psuccess);
double powell_golden_min(int n, int i, double ** u, double ** P,
        double func(double *x,void *),void *,int itmax, double eps);
double ** powell_bracket(int n,int i, double ** u,double **P,
        double func(double *x,void *),void *,int itmax, double eps);

#endif
