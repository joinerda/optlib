#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <optlib.h>

#ifndef HUGE
#define HUGE 1.0e100
#endif

typedef struct {
    int n;
    int fcount;
    double wallcharge;
    double pointcharge;
} Data;

double func(double * x, void * func_data) {
    Data * theData = (Data *)func_data;
    theData->fcount++;
    double retval;
    double wall_charge=theData->wallcharge;
    double point_charge=theData->pointcharge;
    int npoints=theData->n;
    int i,j;
    double xi,yi,xj,yj,dx,dy,dr2,dr;
    double dt,db,dl;

    retval=0.0;
    for(i=0;i<npoints;i++) {
        xi = x[i*2];
        yi = x[i*2+1];
        // out of box penalty
        if(xi<=0.0) return 1.0e50;
        if(xi>=1.0) return 1.0e50;
        if(yi<=0.0) return 1.0e50;
        if(yi>=1.0) return 1.0e50;
        // other charges
        for(j=i+1;j<npoints;j++) {
            xj = x[j*2];
            yj = x[j*2+1];
            dx = xj-xi;
            dy = yj-yi;
            dr2 = dx*dx+dy*dy;
            dr = sqrt(dr2);
            retval += point_charge/dr;
        }
        dt = 1-yi;
        db = yi;
        dl = xi;
        dr = 1-xi;
        retval +=wall_charge*log((sqrt(dt*dt+dr*dr)+dr)/(sqrt(dt*dt+dl*dl)-dl));
        retval +=wall_charge*log((sqrt(db*db+dr*dr)+dr)/(sqrt(db*db+dl*dl)-dl));
        retval +=wall_charge*log((sqrt(dl*dl+db*db)+db)/(sqrt(dl*dl+dt*dt)-dt));
        retval +=wall_charge*log((sqrt(dr*dr+db*db)+db)/(sqrt(dr*dr+dt*dt)-dt));
    }
    return retval;
}


int main(int argc, char ** argv) {
    Data theData;
    int n=10;
    double *guess;
    int i=0;
    int j;
    int seed_offset=0;
    OptLibOpts * theOpts;
    int method=0;

    theData.wallcharge=10.0;
    theData.pointcharge=1.0;
    theData.fcount=0.0;
    if(argc>1) sscanf(argv[1],"%d",&method);
    if(argc>2) sscanf(argv[2],"%d",&n);
    if(argc>3) sscanf(argv[3],"%lf",&theData.wallcharge);
    if(argc>4) sscanf(argv[4],"%lf",&theData.pointcharge);
    theData.n=n/2;
    seed_by_time(0);

    guess = (double *)malloc(sizeof(double)*n);
    
    theOpts = OPTLIB_CreateOpts();
    // initialize guess
    for(i=0;i<n;i++) guess[i] = drand(0.01,0.99);

    // run OPTLIB_Minimize with defaults
    theOpts->method=method;
    theOpts->ga_itmax=1000;
    theOpts->ga_itmin=10;
    theOpts->ga_dominancefactor=0.4;
    theOpts->sa_stepmax=0.75;
    theOpts->sa_nannealers=32;
    OPTLIB_Minimize(n,guess,&func,&theData,theOpts);

    // output
    //for(i=0;i<n;i++) printf("%lf ",guess[i]);
    printf("%lf %d \n",func(guess,&theData),theData.fcount);
    OPTLIB_DestroyOpts(theOpts);

    free(guess);
   
    return 0;
}
