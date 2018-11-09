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

#include "ol_rand_tools.h"

void seed_by_time(int offset) {
    time_t the_time;
    time(&the_time);
    srand((int)the_time+offset);
}

double drand_norm(double xbar, double sigma, double alpha) {
    double omalpha = 1.0-alpha;
    return xbar + sigma*inverf(drand(-omalpha,omalpha));
}

double inverf(double y) {
    // inverse error function computed using Newton's method
    double f,fprime,x;

    x = 0;
    f = erf(x)-y;
    while(fabs(f)>0.0001) {
        fprime = 1.0/atan(1.0)*exp(-x*x);
        x -= f/fprime;
        f = erf(x)-y;
    }

    return x;
}

int irand(int min, int max) {
    return min+(int)((double)rand()/(double)RAND_MAX*(double)(max-min));
}

double drand(double min,double max) {
    return min+((double)rand()/(double)RAND_MAX)*(max-min);
}

void random_direction_subset(int n, int m, double *x) {
    int * deleted;
    int n_deleted=0;
    int direction,exists,i;
    double sum;

    deleted = alloc_iarray(n);
    random_direction(n,x);
    //eliminate n-m directions
    while(n_deleted<n-m) {
        direction = irand(0,n-1);
        exists=0;
        for(i=0;i<n_deleted&&!exists;i++) {
            if(deleted[i]==direction) exists=1;
        }
        if(!exists) {
            deleted[n_deleted]=direction;
            n_deleted++;
            x[direction]=0.0;
        }
    }
    //renormalize
    sum=0.0;
    for(i=0;i<n;i++) {
        sum+=x[i]*x[i];
    }
    sum=1.0/sqrt(sum);
    for (i=0;i<n;i++) {
        x[i]*=sum;
    }
    
    free(deleted);
}

void random_direction(int n, double * x) {
    double length;
    int i;
    do {
        for (i=0;i<n;i++) {
            x[i]=drand(-1,1);
        }
        length = length_darray(n,x);
        if(length!=0.0) {
            for (i=0;i<n;i++) {
                x[i]/=length;
            }
        }
    } while (length==0.0);
}

double drand_updown(double x, double scale) {
    if(drand(0,1)<0.5) {
        return x / drand(0,scale);
    } else {
        return x * drand(0,scale);
    }
}

