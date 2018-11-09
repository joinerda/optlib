#include <optlib.h>
#include <stdlib.h>
#include <stdio.h>
#include "func.h"


int main(int argc, char ** argv) {
    OptLibOpts * theOpts;
    int n=2;
    double  *guess;
    int i;
    func_data_struct func_data;
    func_data.fcount=0;

    printf("OPTLIB DRIVER version %s.%s\n",optlib_VERSION_MAJOR,optlib_VERSION_MINOR);

    guess = (double *)malloc(n*sizeof(double));
 
    seed_by_time(0);

    theOpts = OPTLIB_CreateOpts();
    theOpts->sa_T=5.0;
    theOpts->sa_step=0.5;
    theOpts->sa_stepmax=10.0;
    theOpts->sa_epsilon=1.0e-3;
    theOpts->sa_coolingfactor=0.97;
    theOpts->sa_nannealers=16;
    theOpts->sa_itmax=1000000;
    theOpts->sa_restart=theOpts->sa_itmax+1;
    for (i=0;i<4;i++) {
        func_data.fcount=0;
        func_data.method=i;
        theOpts->method = i;

        guess[0]=drand(0,1);
        guess[1]=drand(0,1);
        OPTLIB_Minimize(n,guess,&func,(void *)&func_data,theOpts);

        if (theOpts->method==OPTLIB_METHOD_GA) printf("METHOD GA\n");
        else if (theOpts->method==OPTLIB_METHOD_SA) printf("METHOD SA\n");
        else if (theOpts->method==OPTLIB_METHOD_MCP) printf("METHOD MCP\n");
        else if (theOpts->method==OPTLIB_METHOD_POWELL)
            printf("METHOD POWELL\n");
        printf("SOLUTION %g %g\n",guess[0], guess[1]);
        printf("TOTAL FUNCTION CALLS %d\n",func_data.fcount);
    }

    OPTLIB_DestroyOpts(theOpts);
    free(guess);
}
