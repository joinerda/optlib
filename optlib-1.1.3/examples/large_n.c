#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <optlib.h>

typedef struct {
    int n;
    int m;
    int count;
    double gamma;
} Data;

// determination of vertical miss distance squared,
//  assuming target is 100 m away and 50 m up.
double func(double *x, void * func_data) {
    Data * theData = (Data *)func_data;
    int n = theData->n;
    int m = theData->m;
    theData->count++;
    double gamma = theData->gamma;
    double r;
    double f;
    int i;

    r = 0.0;
    for(i=0;i<m;i++) {
        r+=(x[i]-0.5)*(x[i]-0.5);
    }
    r = sqrt(r);

    f = cos((double)n*M_PI*r)*exp(-gamma*r*r);
    return r;

}

int main(int argc, char ** argv) {
    int m=8;
    int n=9;
    double gamma=0.15;
    double *guess;
    int i;
    int j;
    int seed_offset=0;
    OptLibOpts * theOpts;
    Data theData;


    if(argc>1) sscanf(argv[1],"%d",&m);
    theData.n=n;
    theData.m=m;
    theData.gamma=gamma;

    guess = (double *)malloc(sizeof(double)*m);

    // always seed your stochastic models
    if(argc>1) sscanf(argv[1],"%d",&seed_offset);
    seed_by_time(seed_offset);
    
    theOpts = OPTLIB_CreateOpts();
    //for(j=0;j<m;j++) printf("x%d ",j);
    printf("e count\n");
    for (i=0;i<10;i++) {
        theData.count=0;
        // initialize guess
        for(j=0;j<m;j++) guess[j]=drand(0.0,1.0);
    
        // run OPTLIB_Minimize with some defaults overridden
        theOpts->ga_itmax=100000;
        theOpts->ga_dominancefactor=0.4;
        //theOpts->method=OPTLIB_METHOD_SA;
        theOpts->sa_coolingfactor=0.96;
        theOpts->sa_stepmax=0.75;
        OPTLIB_Minimize(m,guess,&func,&theData,theOpts);
    
        // output
        //for(j=0;j<m;j++) printf("%lf ",guess[j]);
        printf("%5.3lf %d\n",func(guess,&theData), theData.count);
    }
    OPTLIB_DestroyOpts(theOpts);

   
    free(guess);
    return 0;
}
