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
#include "func.h"

double func(double * par, void * func_data_void) {
    func_data_struct * func_data;
    int method;
    //double x = (atan(par[0])+M_PI/2)/M_PI;
    //double y = (atan(par[1])+M_PI/2)/M_PI;
    double x;
    double y;
    double retval;
    double rr;
    double sigma2=.15;
    int n=9;
    double penalty=0.0;

    x = par[0];
    y = par[1];

    func_data = (func_data_struct *)func_data_void;
    func_data->fcount++;
    method = func_data->method;

/*
    retval = 16.0*x*y*(1.0-x)*(1.0-y)*sin(n*M_PI*x)*sin(n*M_PI*y);
    retval *= retval;
*/

    penalty=0.0;
    if(method==1) {
        if(x<0.0||x>1.0) penalty += pow((5.0+fabs(x)),2.0);
        if(y<0.0||y>1.0) penalty += pow((5.0+fabs(x)),2.0);
    }

// pikaia test example
    rr = sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5));
    retval = pow(cos(rr*n*M_PI),2.0)*exp(-rr*rr/sigma2);

    return -retval+penalty;

}

