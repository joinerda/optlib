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


#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "ol_mem.h"
#include "ol_fw_mpi.h"
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

#ifdef USE_OPTMPI
void fw_worker_main(fw_struct * func_struct) {
    int done=FALSE;
    int mesg;
    double *x;
    double y;
    MPI_Status status;

    x = alloc_darray(func_struct->n);

    while(!done) {
        MPI_Recv(&mesg,1,MPI_INT,0,FW_MESSAGE,func_struct->comm,&status);
        if (mesg==FW_MSG_WORK) {
            MPI_Recv(x,func_struct->n,
                MPI_DOUBLE,0,FW_DATA_IN,func_struct->comm,&status);
            y = func_struct->func(x,func_struct->func_data);
            MPI_Send(&y,1,MPI_DOUBLE,0,FW_DATA_OUT,func_struct->comm);
        } else {
            done=TRUE;
        }
    }

    free(x);
}

int fw_next_worker(int worker, int size, int foreman_works) {
    if(foreman_works) {
        return (worker+1)%size;
    } else {
        return (worker%(size-1))+1;
    }
}

void fw_send_job(fw_struct * func_struct,
        double *x,int *worker, MPI_Request * prequest) {
    int mesg=FW_MSG_WORK;
    if(*worker>0&&func_struct->size>1) {
        MPI_Send(&mesg,1,MPI_INT,*worker,FW_MSG_WORK,func_struct->comm);
        MPI_Isend(x,func_struct->n,MPI_DOUBLE,
            *worker,FW_DATA_IN,func_struct->comm,prequest);
    }
    *worker = fw_next_worker(*worker,func_struct->size,
        func_struct->foreman_works);
}

void fw_shutdown_worker(int worker,MPI_Comm comm) {
    int mesg=FW_MSG_QUIT;
    MPI_Send(&mesg,1,MPI_INT,worker,FW_MESSAGE,comm);
}

void fw_shutdown_workers(int size,MPI_Comm comm) {
    int worker;
    for(worker=1;worker<size;worker++) {
        fw_shutdown_worker(worker,comm);
    }
}

double fw_get_job(fw_struct * func_struct,
        double *x, MPI_Request * prequest, int *worker) {
    MPI_Status status;
    double y;
    if(*worker>0&&func_struct->size>1) {
        MPI_Wait(prequest,&status);
        MPI_Recv(&y,1,MPI_DOUBLE,
            *worker,FW_DATA_OUT,func_struct->comm,&status);
    } else {
        y = func_struct->func(x,func_struct->func_data);
    }
    *worker = fw_next_worker(*worker,func_struct->size,
        func_struct->foreman_works);
    return y;
}

void fw_simd_func(double **x,double *y,int n_jobs,fw_struct * func_struct) {
    int worker,i;
    MPI_Request * request;

    request = (MPI_Request*)malloc(sizeof(MPI_Request)*n_jobs);

    // prep all workers and send jobs
    worker=1;
    for(i=0;i<n_jobs;i++) {
        fw_send_job(func_struct,x[i],&worker,&(request[i]));
    }

    // wait for workers to finish
    worker=1;
    for(i=0;i<n_jobs;i++) { 
        y[i] = fw_get_job(func_struct,x[i],&(request[i]),&worker);
    }

    free(request);
    return;
}
#else
// placeholder for serial version of fw_mpi codes
void fw_simd_func(double **x,double *y,int n_jobs,fw_struct * func_struct) {
    int i;

    for(i=0;i<n_jobs;i++) {
        y[i]=func_struct->func(x[i],func_struct->func_data);
    }

    return;
}
#endif

