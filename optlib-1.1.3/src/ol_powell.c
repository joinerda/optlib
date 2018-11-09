/**********************************************
*  Copyright 2007-2008, David Joiner
*
*  Powell's method data fit
*
**********************************************/



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
#include "ol_mem.h"
#include "ol_powell.h"
#include "ol_rand_tools.h"

#define GOLDC 0.38197
#define GOLDR (1.0-GOLDC)
#define GOLD 1.618034
#define TRYNEWTON 0
#define EXPANSION_SCALE 0.1

/* Copyright 2008 David Joiner */

/* Algorithms taken from Numerical Recipes for powell method
   and for golden rule method of finding linear min */

double powell_directional_deriv(
        int n, int i, double ** u,double * x,double func(double * x, void *),
        void * func_data) {

    double *step;
    double *x_plus;
    double *x_minus;
    double h;
    int l;
    double return_value;

    step=alloc_darray(n);
    x_plus=alloc_darray(n);
    x_minus=alloc_darray(n);

    h = 0.0;
    for(l=0;l<n;l++) {
        step[l] = 0.001*u[l][i];
        h += step[l]*step[l];
        x_plus[l] = x[l]+step[l];
        x_minus[l] = x[l]-step[l];
    }
    h = 2.0*sqrt(h);

    return_value = (func(x_plus,func_data)-func(x_minus,func_data))/h;

    free_darray(step);
    free_darray(x_plus);
    free_darray(x_minus);

    return return_value;
}

double powell_directional_deriv2(
        int n, int i, double ** u,double * x,double func(double * x,void *),
        void * func_data) {

    double *step;
    double *x_plus;
    double *x_minus;
    double h;
    int l;
    double return_value;

    step=alloc_darray(n);
    x_plus=alloc_darray(n);
    x_minus=alloc_darray(n);

    h = 0.0;
    for(l=0;l<n;l++) {
        step[l] = 0.001*u[l][i];
        h += step[l]*step[l];
        x_plus[l] = x[l]+step[l];
        x_minus[l] = x[l]-step[l];
    }
    h = 2.0*sqrt(h);

    return_value =
        (powell_directional_deriv(n,i,u,x_plus,func,func_data)-
        powell_directional_deriv(n,i,u,x_minus,func,func_data))/h;

    free_darray(step);
    free_darray(x_plus);
    free_darray(x_minus);

    return return_value;
}

double vector_distance(int n, double * x, double * y) {
    double length;
    int i;

    length = 0.0;
    for(i=0;i<n;i++) {
        length += pow(x[i]-y[i],2.0);
    }
    length = sqrt(length);
    return length;
}

double vector_length(int n, double * x) {
    double length;
    int i;

    length = 0.0;
    for(i=0;i<n;i++) {
        length += pow(x[i],2.0);
    }
    length = sqrt(length);
    return length;
}

double powell_golden_min(int n, int i, double ** u, double ** P,
        double func(double *x,void *),void * func_data,
        int itmax, double eps) {

    double ** bracket;
    double * x0;
    double * x1;
    double * x2;
    double * x3;
    double * x_old;
    double f_old;
    double f0,f1,f2,f3;
    double left,right;
    double delta_f;
    int l;

    x0 = alloc_darray(n);
    x1 = alloc_darray(n);
    x2 = alloc_darray(n);
    x3 = alloc_darray(n);
    x_old = alloc_darray(n);

    if(i==0) {
        for(l=0;l<n;l++) {
            x_old[l] = P[l][n];
        }
    } else {
        for(l=0;l<n;l++) {
            x_old[l] = P[l][i-1];
        }
    }
    f_old = func(x_old,func_data);
    bracket = powell_bracket(n,i,u,P,func,func_data,itmax,eps);

    left = 0.0;
    right = 0.0;
    for (l=0;l<n;l++) {
        x0[l] = bracket[l][0];
        x3[l] = bracket[l][2];
        left += pow(bracket[l][1]-x0[l],2.0);
        right += pow(x3[l]-bracket[l][1],2.0);
    }
    f0 = bracket[n][0];
    f3 = bracket[n][2];
    left = sqrt(left);
    right = sqrt(right);

    // determine larger interval
    if(left>right) {
        for(l=0;l<n;l++) {
            x1[l]=bracket[l][1];
            x2[l]=x1[l]+GOLDC*(x3[l]-x1[l]);
        }
    } else {
        for(l=0;l<n;l++) {
            x2[l]=bracket[l][1];
            x1[l]=x2[l]-GOLDC*(x2[l]-x0[l]);
        }
    }
    f1=func(x1,func_data);
    f2=func(x2,func_data);
    // check each component individually in stopping condition
    while(vector_distance(n,x3,x0) >
            eps*(vector_length(n,x1)+vector_length(n,x2))) {
        if(f2<f1) {
            for(l=0;l<n;l++) {
                x0[l]=x1[l];
                x1[l]=x2[l];
                x2[l]=GOLDR*x1[l]+GOLDC*x3[l];
            }
            f1=f2;
            f2=func(x2,func_data);
        } else {
            for(l=0;l<n;l++) {
                x3[l]=x2[l];
                x2[l]=x1[l];
                x1[l]=GOLDR*x2[l]+GOLDC*x0[l];
            }
            f2=f1;
            f1=func(x1,func_data);
        }
    }
    if(f1<f2) {
        // set P[i]
        // determine change
        delta_f = f_old-f1;
        for(l=0;l<n;l++) {
            P[l][i] = x1[l];
        }
    } else {
        // set P[i]
        // determine change
        delta_f = f_old-f2;
        for(l=0;l<n;l++) {
            P[l][i] = x2[l];
        }
    }
    delta_f = fabs(delta_f);
    
    free_dmatrix(bracket);
    free(x0);
    free(x1);
    free(x2);
    free(x3);
    free(x_old);

    //if(delta_f<0.0) printf("WARNING! %10.3e %10.3e %10.3e \n",f_old,f1,f2);

    return delta_f;
}

void print_vector(const char * mesg,int n,double * x) {
    int l;
    printf("%s\t",mesg);
    for (l=0;l<n;l++) printf("%10.3e\t",x[l]);
    printf("\n");
}

double ** powell_bracket(int n,int i, double ** u,double **P,
        double func(double *x,void *),void * func_data,int itmax, double eps) {
    double *x_left;
    double *x_right;
    double *x_mid;
    double ** bracket;
    double func_mid,func_left,func_right;
    int l;
    int im1;
    double factor;
    int iter;
    double deriv;

    x_left=alloc_darray(n);
    x_right=alloc_darray(n);
    x_mid=alloc_darray(n);
    bracket = alloc_dmatrix(n+1,3);

    if (i>0) {
        im1=i-1;
        for (l=0;l<n;l++) {
            x_mid[l]=P[l][im1];
        }
    } else {
        im1=n-1;
        for (l=0;l<n;l++) {
            x_mid[l]=P[l][n];
        }
    }
    factor = GOLD;

    // start with current guess. Expand outwards for two outer
    // guesses until 
    // use directional derivative to decide which way to expand.
    deriv = powell_directional_deriv(n,im1,u,x_mid,func,func_data);
    for (l=0;l<n;l++) {
        if (deriv<=eps) {
            x_left[l]=x_mid[l]-vector_length(n,x_mid)*u[l][im1]*GOLD*EXPANSION_SCALE;
            x_right[l]=x_mid[l]+vector_length(n,x_mid)*u[l][im1]*GOLDC*EXPANSION_SCALE;
        } 
        if (deriv>=-eps) {
            x_right[l]=x_mid[l]+vector_length(n,x_mid)*u[l][im1]*GOLD*EXPANSION_SCALE;
            x_left[l]=x_mid[l]-vector_length(n,x_mid)*u[l][im1]*GOLDC*EXPANSION_SCALE;
        }
    }

    // test to see if minimum has been bracketed (i.e. is the
    //   function on both sides greater than the middle?
    func_mid = func(x_mid,func_data);
    func_left = func(x_left,func_data);
    func_right = func(x_right,func_data);
    iter=0;
    while(((func_left<func_mid)||(func_right<func_mid)) && iter<itmax) {
        // determine downhill side
        if (func_left<func_right) {
            // increase left ratio
            for(l=0;l<n;l++) {
                x_left[l] = x_left[l] - 0.5*(x_right[l]-x_left[l]);
            }
            func_left = func(x_left,func_data);
        } else {
            // increase right ratio
            for(l=0;l<n;l++) {
                x_right[l] = x_right[l] + 0.5*(x_right[l]-x_left[l]);
            }
            func_right = func(x_right,func_data);
        }
        // update mid point
        if (func_left<func_right) {
            for(l=0;l<n;l++) {
                x_mid[l] = x_left[l] + GOLDC*(x_right[l]-x_left[l]);
            }
        } else {
            for(l=0;l<n;l++) {
                x_mid[l] = x_right[l] - GOLDC*(x_right[l]-x_left[l]);
            }
        }
        func_mid = func(x_mid,func_data);
        iter++;
    }
    if(iter>=itmax*10) {
        printf("WARNING: itmax exceeded in powell_bracket\n");
    }
    for(l=0;l<n;l++) {
        bracket[l][0] = x_left[l];
        bracket[l][1] = x_mid[l];
        bracket[l][2] = x_right[l];
    }
    bracket[n][0] = func_left;
    bracket[n][1] = func_mid;
    bracket[n][2] = func_right;

    free(x_left);
    free(x_right);
    free(x_mid);
    return bracket;
}

// update P[i] along u[i]
double powell_findmin(int n,int i,double ** u,
        double ** P,double func(double * x,void *),
        void * func_data, int itmax, double eps) {
    int isuccess;
    double test_newton;
    if(TRYNEWTON) {
        test_newton=powell_newton_min(n,i,u,P,func,func_data,
            itmax,eps,&isuccess);
        if(isuccess) {
            return test_newton;
        } else {
            return powell_golden_min(n,i,u,P,func,func_data,
                itmax,eps);
        }
    } else {
        return powell_golden_min(n,i,u,P,func,func_data,
            itmax,eps);
    }
}

double powell_newton_min(int n,int i,double ** u,
        double ** P,double func(double * x,void *),
        void * func_data, int itmax, double eps,
        int * psuccess) {
    // Use Newton's Method
    //  calculate first and second derivative, iterate x
    //  until first derivative equals zero.
    double *x_new;
    double *x_old;
    double *step;
    double deriv,deriv2;
    double h;
    int l,iter,done;
    int im1;
    double init_f,final_f;

    x_new=alloc_darray(n);
    x_old=alloc_darray(n);
    step=alloc_darray(n);

    if (i>0) {
        im1=i-1;
        for (l=0;l<n;l++) {
            x_old[l]=P[l][im1];
        }
    } else {
        im1=n-1;
        for (l=0;l<n;l++) {
            x_old[l]=P[l][n];
        }
    }
    init_f = func(x_old,func_data);
    
    done=0;
    iter=0;
    *psuccess=0;
    deriv=1.0;
    deriv2=0.0;
    while(!done&&iter<itmax) {
        deriv = powell_directional_deriv(n,im1,u,x_old,
            func,func_data);
        deriv2 = powell_directional_deriv2(n,im1,u,x_old,
            func,func_data);
        for (l=0;l<n;l++) {
            step[l] = -u[l][im1]*deriv/deriv2;
            x_new[l] = x_old[l] + step[l];
            x_old[l] = x_new[l];
        }
        h=0.0;
        for (l=0;l<n;l++) h += step[l]*step[l];
        h = sqrt(h);
        if (h<eps) done=1;
        iter++;
    }
    if(iter>=itmax) {
        printf("WARNING: # iterations exceeded in findmin\n");
    }

    for (l=0;l<n;l++) P[l][i]=x_new[l];
    final_f = func(x_old,func_data);

    free_darray(x_new);
    free_darray(x_old);
    free_darray(step);

    if(deriv2>0.0) *psuccess=1;

    return init_f-final_f;
}

#ifdef USE_OPTMPI
struct { 
    double val; 
    int rank; 
} in_minloc, out_minloc;
#endif

void mc_powell(int n_trials, double scale,
        double * x_attempts, double * x_values, double * f_values,
        int n,double *x,double func(double *,void *),
        void * func_data,void output(double *, void *, const char *),
        int itmax,double eps) {
    double * return_value = NULL;
    double fvalue_min = 0.0;
    double * xvalue_min;
    double * xorig;
#ifdef USE_OPTMPI
    double * buffer;
#endif
    double ftest;
    int i,j,ij;
    int size=1;

#ifdef USE_OPTMPI
    buffer = alloc_darray(n*n_trials);
#endif
    xvalue_min = alloc_darray(n);
    xorig = alloc_darray(n);

    for(i=0;i<n;i++) xorig[i]=x[i];
    if(x_attempts!=NULL) for(ij=0;ij<n_trials*n;ij++) x_attempts[ij]=0.0;
    if(x_values!=NULL) for(ij=0;ij<n_trials*n;ij++) x_values[ij]=0.0;
    if(f_values!=NULL) for(j=0;j<n_trials;j++)f_values[j]=0.0;

    // use initial guess to create minimum value
    for(i=0;i<n;i++) {
        xvalue_min[i]=x[i];
        if(x_attempts!=NULL) x_attempts[0*n+i]=x[i];
    }
    powell(n,xvalue_min,func,func_data,output,itmax,eps);
    for(i=0;i<n;i++) {
        if(x_values!=NULL) x_values[0*n+i]=xvalue_min[i];
    }
    fvalue_min = func(xvalue_min,func_data);
    if(f_values!=NULL) f_values[0]=fvalue_min;
    
    if(scale==0.0) scale=10.0;

#ifdef USE_OPTMPI
    MPI_Comm_size(mcp_mpi_comm,&size);
#endif
    
    for(j=1;j<n_trials;j+=size) {
        for(i=0;i<n;i++) {
            if(xorig[i] != 0.0) {
                x[i] = xorig[i] * (1.0 + scale *drand(-1.0,1.0));
            } else {
                x[i] = xorig[i] + scale * drand(-1.0,1.0);
            }
        }
        for(i=0;i<n;i++) {
            if(x_attempts!=NULL) x_attempts[j*n+i]=x[i];
        }
        powell(n,x,func,func_data,output,itmax,eps);
        ftest = func(x,func_data);
        for(i=0;i<n;i++) {
            if(x_values!=NULL) x_values[j*n+i]=x[i];
        }
        if(f_values!=NULL) f_values[j]=ftest;
        if(ftest<fvalue_min) {
            fvalue_min = ftest;
            for(i=0;i<n;i++) xvalue_min[i] = x[i];
        }
    }
#ifdef USE_OPTMPI
    if(x_attempts!=NULL) {
        MPI_Allreduce(x_attempts,buffer,n*n_trials,MPI_DOUBLE,MPI_SUM,
            mcp_mpi_comm);
        for(ij=0;ij<n*n_trials;ij++) x_attempts[ij]=buffer[ij];
    }
    if(x_values!=NULL) {
        MPI_Allreduce(x_values,buffer,n*n_trials,MPI_DOUBLE,MPI_SUM,
            mcp_mpi_comm);
        for(ij=0;ij<n*n_trials;ij++) x_values[ij]=buffer[ij];
    }
    if(f_values!=NULL) {
        MPI_Allreduce(f_values,buffer,n_trials,MPI_DOUBLE,MPI_SUM,
            mcp_mpi_comm);
        for(j=0;j<n_trials;j++) f_values[j]=buffer[j];
    }
    in_minloc.val = fvalue_min;
    in_minloc.rank = mcp_mpi_rank;
    MPI_Allreduce(&in_minloc,&out_minloc,1,
        MPI_DOUBLE_INT,MPI_MINLOC,mcp_mpi_comm);
    MPI_Bcast(&fvalue_min,1,MPI_DOUBLE,out_minloc.rank,mcp_mpi_comm);
    MPI_Bcast(xvalue_min,n,MPI_DOUBLE,out_minloc.rank,mcp_mpi_comm);
#endif

    for(i=0;i<n;i++) x[i]=xvalue_min[i];

    free(xvalue_min);
    free(xorig);
#ifdef USE_OPTMPI
    free(buffer);
#endif

    return;
        
}

void powell(int n,double *x,double func(double *,void *),
        void * func_data,void output(double *, void *, const char *),
        int itmax,double eps) {

    double **P;
    double **u;
    int i,j,l;
    int iter;
    double h;
    double *delta_f;
    double *pzero;
    double *pn;
    double *twopn_minus_pzero;
    double fzero,fn,fe;
    double delta_max;
    double i_max;
    int done;

    P = alloc_dmatrix(n,n+1);
    u = alloc_dmatrix(n,n);
    delta_f = alloc_darray(n);
    pzero = alloc_darray(n);
    pn = alloc_darray(n);
    twopn_minus_pzero = alloc_darray(n);

    // Set Pzero to the initial guess
    for(l=0;l<n;l++) P[l][0] = x[l];

    // Set initial directions to basis vectors
    for (l=0;l<n;l++) {
        for(j=0;j<n;j++) {
            if(l!=j) {
                u[l][j]=0.1;
            } else {
                u[l][j]=0.0;
            }
        }
    }
    
    done=0;
    iter=0;
    while(!done && iter<itmax) {
        // move points along each directions to find minimum.
        delta_max=-1.0;
        i_max=-1;
        for (i=1;i<n+1;i++) {
            delta_f[i-1] = powell_findmin(n,i,u,P,func,
                func_data,itmax,eps);
            if(delta_f[i-1]>delta_max) {
                delta_max=delta_f[i-1];
                i_max=i-1;
            }
        }
        for(l=0;l<n;l++) {
            pzero[l]=P[l][0];
            pn[l]=P[l][n];
            twopn_minus_pzero[l]=2.0*pn[l]-pzero[l];
        }
        fzero = func(pzero,func_data);
        if(POWELL_VERBOSE) {
            if(output!=NULL) {
                output(pzero,func_data,"POWELL X ");
            } else {
                print_darray(n,pzero,"POWELL X ");
            }
            printf("LOW ENERGY  = %10.3e\n",fzero);
        }
        fn = func(pn,func_data);
        fe = func(twopn_minus_pzero,func_data);

        if(!(fe>=fzero)||
                !((2.0*(fzero-2.0*fn+fe)*pow((fzero-fn)-delta_max,2.0))>=
                  pow(fzero-fe,2.0)*delta_max)) {
            // reset directions, discard last best direction
            for (i=i_max;i<n-1;i++) {
                for(l=0;l<n;l++) {
                    u[l][i]=u[l][i+1];
                }
            }
            h=0.0;
            for (l=0;l<n;l++) {
                u[l][n-1]=P[l][n]-P[l][0];
                h+=u[l][n-1]*u[l][n-1];
            }
            // normalize u
            h = sqrt(h);
            for (l=0;l<n;l++) {
                u[l][n-1]=u[l][n-1]/h;
            }
        }
        powell_findmin(n,0,u,P,func,func_data,itmax,eps);
        h = 0.0;
        for (l=0;l<n;l++) {
            h += pow(P[l][0]-pzero[l],2.0);
        }
        h = sqrt(h);
        if (h<eps) {
            done=1;
        }
        iter++;
    }
    for(l=0;l<n;l++) x[l]=P[l][0];


    free_dmatrix(P);
    free_dmatrix(u);
    free_darray(delta_f);
    free_darray(pzero);
    free_darray(pn);
    free_darray(twopn_minus_pzero);

}

