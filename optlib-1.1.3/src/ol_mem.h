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

#ifndef MEM_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#undef min
#undef max
#define min(X,Y) ((X) < (Y) ? (X) : (Y))
#define max(X,Y) ((X) > (Y) ? (X) : (Y))

double length_darray(int n,double *x);
double length_farray(int n,float *x);
double * alloc_darray(int n);
double ** alloc_dmatrix(int n,int m);
void free_dmatrix(double ** dmatrix);
void free_darray(double * darray);
float * alloc_farray(int n);
float ** alloc_fmatrix(int n,int m);
void free_fmatrix(float ** fmatrix);
void free_farray(float * farray);
int * alloc_iarray(int n);
int ** alloc_imatrix(int n,int m);
void free_imatrix(int ** imatrix);
void free_iarray(int * iarray);
void print_darray(int n, double * x, const char * string);
void print_darrayt(int n, double * x, const char * string);

#define MEM_H
#endif

