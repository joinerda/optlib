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

#include "ol_mem.h"

double length_darray(int n, double *x) {
    int i;
    double sum;

    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    return sqrt(sum);
}

double length_farray(int n, float *x) {
    int i;
    float sum;

    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    return sqrt(sum);
}

double * alloc_darray(int n) {
    int i;
    double *darray = (double *) malloc(sizeof(double)*n);
    if(darray==NULL) {
        printf("WARNING: failed to allocate memory in alloc_array!\n");
    }
    for (i=0;i<n;i++) {
        darray[i]=0.0;
    }
    return darray;
}

double ** alloc_dmatrix(int n,int m) {
    int i,j;
    double ** dmatrix;
    double * large_array;
    large_array=alloc_darray(n*m);
    dmatrix = (double **) malloc(sizeof(double*)*n); 
    dmatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        //dmatrix[i] = &large_array[i*m];
        dmatrix[i] = dmatrix[i-1]+m;
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            dmatrix[i][j]=0.0;
        }
    }
    return dmatrix;
}

void free_dmatrix(double ** dmatrix) {
    free(*dmatrix);
    free(dmatrix);
}

void free_darray(double * darray) {
    if(darray!=NULL) free(darray);
}

float * alloc_farray(int n) {
    int i;
    float *farray = (float *) malloc(sizeof(float)*n);
    for (i=0;i<n;i++) {
        farray[i]=0.0;
    }
    return farray;
}

float ** alloc_fmatrix(int n,int m) {
    int i,j;
    float ** fmatrix;
    float * large_array;
    large_array=alloc_farray(n*m);
    fmatrix = (float **) malloc(sizeof(float*)*n); 
    fmatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        fmatrix[i] = &large_array[i*m];
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            fmatrix[i][j]=0.0;
        }
    }
    return fmatrix;
}

void free_fmatrix(float ** fmatrix) {
    free(*fmatrix);
    free(fmatrix);
}

void free_farray(float * farray) {
    free(farray);
}

int * alloc_iarray(int n) {
    int i;
    int *iarray = (int *) malloc(sizeof(int)*n);
    for (i=0;i<n;i++) {
        iarray[i]=0.0;
    }
    return iarray;
}

int ** alloc_imatrix(int n,int m) {
    int i,j;
    int ** imatrix;
    int * large_array;
    large_array=alloc_iarray(n*m);
    imatrix = (int **) malloc(sizeof(int*)*n); 
    imatrix[0] = &large_array[0];
    for(i=1;i<n;i++) {
        imatrix[i] = &large_array[i*m];
    }
    for(i=0;i<n;i++) {
        for(j=0;j<m;j++) {
            imatrix[i][j]=0.0;
        }
    }
    return imatrix;
}

void free_imatrix(int ** imatrix) {
    free(*imatrix);
    free(imatrix);
}

void free_iarray(int * iarray) {
    free(iarray);
}

void print_darray(int n, double * x, const char * string) {
    int i;
    printf("%s",string);
    for(i=0;i<n;i++) printf("[%d] %10.3e\t",i,x[i]);
    printf("\n");

}
void print_darrayt(int n, double * x, const char * string) {
    int i;
    printf("%s",string);
    for(i=0;i<n;i++) printf("TRANSFORMED [%d] %10.3e\t",i,(atan(x[i])+M_PI/2)/M_PI);
    printf("\n");

}
