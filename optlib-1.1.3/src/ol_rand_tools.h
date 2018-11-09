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

#ifndef RANDTOOLS
#define RANDTOOLS
#include <time.h>
#include <stdlib.h>
#include "ol_mem.h"

void seed_by_time(int offset);
double drand(double min,double max);
double drand_norm(double xbar,double sigma,double alpha);
double inverf(double);
int irand(int min,int max);
void random_direction(int n, double * x);
void random_direction_subset(int n, int m, double * x);
double drand_updown(double x, double scale);


#endif
