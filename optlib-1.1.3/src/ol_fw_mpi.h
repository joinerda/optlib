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

#ifndef FW_MPI
#define FW_MPI

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ol_mem.h"
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

#define FW_MESSAGE 1
#define FW_DATA_IN 2
#define FW_DATA_OUT 3

#define FW_MSG_WORK 1
#define FW_MSG_QUIT 2

#ifndef FALSE
#define FALSE 0
#endif

#ifndef TRUE
#define TRUE 1
#endif

typedef struct {
    double (*func)(double *,void *);
    void *func_data;
    int size;
    int rank;
#ifdef USE_OPTMPI
    MPI_Comm comm;
#endif
    int n;
    int foreman_works;
} fw_struct;

#ifdef USE_OPTMPI
void fw_worker_main(fw_struct * func_struct);
int fw_next_worker(int worker, int size, int foreman_works);
void fw_send_job(fw_struct * func_struct,
        double *x,int *worker, MPI_Request * prequest);
void fw_shutdown_worker(int worker,MPI_Comm comm);
void fw_shutdown_workers(int size,MPI_Comm comm);
double fw_get_job(fw_struct * func_struct,
        double *x, MPI_Request * prequest, int *worker);
#endif
void fw_simd_func(double **x,double *y,int n_jobs,
        fw_struct * func_struct);

#endif
