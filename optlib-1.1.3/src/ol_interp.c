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
#include "ol_interp.h"

double log_interp(double x,int n,double * x_data,double * y_data) {
    int i;
    double y;
    double * log_x_data;
    double * log_y_data;

    log_x_data = alloc_darray(n);
    log_y_data = alloc_darray(n);

    for(i=0;i<n;i++) {
        log_x_data[i] = log(x_data[i]);
        log_y_data[i] = log(y_data[i]);
    }
    x = log(x);
    y = lin_interp(x,n,log_x_data,log_y_data);
    y = exp(y);

    free(log_x_data);
    free(log_y_data);

    return y;
}

double lin_interp(double x,int n,double * x_data,double * y_data) {
    // determine if x is inside of x_data
    int i;
    double y;
    if(x_data[0]<x_data[n-1]) {
        // increasing
        if(x<x_data[0]||x>x_data[n-1]) {
            // out of range
            printf("WARNING: OUT OF RANGE IN INTERPOLATION 1\n");
            return 0.0;
        }
    } else {
        // decreasing
        if(x>x_data[0]||x<x_data[n-1]) {
            // out of range
            printf("WARNING: OUT OF RANGE IN INTERPOLATION 2\n");
            printf("WARNING: %10.3e %10.3e %10.3e\n",x_data[0],x,x_data[n-1]);
            return 0.0;
        }
    }
    y=0.0;
    for(i=0;i<n-1;i++) {
        if((x<x_data[i]&&x>=x_data[i+1])||(x>x_data[i]&&x<=x_data[i+1])) {
            // found interval, interpolate
            if(fabs(x-x_data[i])<=fabs(x-x_data[i+1])) {
                y = y_data[i]+(y_data[i+1]-y_data[i])/(x_data[i+1]-x_data[i])*
                    (x-x_data[i]);
            } else {
                y = y_data[i+1]+(y_data[i+1]-y_data[i])/(x_data[i+1]-x_data[i])*
                    (x-x_data[i+1]);
            }
            break;
        }
    }
    return y;
}

double poly_interp(double x,int n,double * x_data,double * y_data) {
    // determine if x is inside of x_data
    int i,j,k;
    double y;
    double sum,product;
    int start,npoints;
    if(x_data[0]<x_data[n-1]) {
        // increasing
        if(x<x_data[0]||x>x_data[n-1]) {
            // out of range
            printf("WARNING: OUT OF RANGE IN INTERPOLATION\n");
            return 0.0;
        }
    } else {
        // decreasing
        if(x>x_data[0]||x<x_data[n-1]) {
            // out of range
            printf("WARNING: OUT OF RANGE IN INTERPOLATION\n");
            return 0.0;
        }
    }
    y=0.0;
    for(i=0;i<n-1;i++) {
        if((x<x_data[i]&&x>=x_data[i+1])||(x>x_data[i]&&x<=x_data[i+1])) {
            // found interval, interpolate
            if(i==0) {
                start=i;
                npoints=3;
            } else if(i==n-2) {
                start=i-1;
                npoints=3;
            } else {
                start=i-1;
                npoints=4;
            }
            sum=0.0;
            for(k=start;k<start+npoints;k++) {
                product=y_data[k];
                for(j=start;j<start+npoints;j++) {
                    if(j!=k) {
                        product *= (x-x_data[j])/(x_data[k]-x_data[j]);
                    }
                }
                sum=sum+product;
            }
            y = sum;
            break;
        }
    }
    return y;
}
