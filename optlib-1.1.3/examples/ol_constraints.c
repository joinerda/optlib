#include <optlib_mpi.h>
#include <stdlib.h>
#include <stdio.h>

double func(double *x, void * func_data) {
    int * fcount = (int *)func_data;
    double y;
    int a;
    int d;
    double z;

    (*fcount)++;

    y = x[0];

    a = (int)x[1]; // value constrained to an integer

    // enumerated values, 
    //   penalize anything out of range, use mod
    //   within range
    if(x[2]<0.0) return 100*(10-x[2])+1000.0;
    if(x[2]>=3.0) return 100*(10.0+x[2])+1000.0;

    d = ((int)x[2])%3; // enum between 0 and 2

    // constraining of positive-only values can easily
    // be done by transforming the variable to be
    // minimized
    z = exp(x[3]);  // x[3] must be positive, so
                    // optimize on log of value instead
                    // of value...inverse of log is exp.

    return (y-0.5)*(y-0.5)+
        (a-5)*(a-5)+
        (d-1)*(d-1)+   
        (sqrt(z)-2)*(sqrt(z)-2);
}


int main(int argc, char ** argv) {
    OptLibOpts * theOpts;
    int n=4;
    double guess[4];
    int i;
    int fcount=0;
    int rank;
    int size;
    int buffer;

    MPI_Init(&argc,&argv);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm_size(MPI_COMM_WORLD,&size);

    if(rank==0) printf("OPTLIB DRIVER version %s.%s\n",optlib_VERSION_MAJOR,optlib_VERSION_MINOR);

    seed_by_time(rank);

    if(rank==0) {
        guess[0]=drand(-10.0,10.0);
        guess[1]=drand(-10.0,10.0);
        guess[2]=drand(0.0,3.0);
        guess[3]=log(drand(0.1,100.0));
    }
    MPI_Bcast(guess,n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    theOpts = OPTLIB_CreateOpts();
    theOpts->method = OPTLIB_METHOD_SA;
    theOpts->sa_nannealers=1;
    theOpts->sa_coolingfactor=0.90;
    OPTLIB_Minimize(n,guess,&func,&fcount,theOpts);
    MPI_Reduce(&fcount,&buffer,1,MPI_INT,MPI_SUM,0,MPI_COMM_WORLD);
    if(rank==0) {
        fcount=buffer;
        printf("TOTAL FUNCTION CALLS %d\n",fcount);
    }

    // cleanup final result with Powell's method
    theOpts->method = OPTLIB_METHOD_POWELL;
    theOpts->powell_epsilon=1.0e-10;
    fcount=0;
    //OPTLIB_Minimize(n,guess,&func,&fcount,theOpts);

    if(rank==0) printf("SOLUTION %g %d %d %g\n",guess[0], (int)guess[1], (int)guess[2],exp(guess[3]));
    //printf("TOTAL FUNCTION CALLS %d\n",fcount);

    OPTLIB_DestroyOpts(theOpts);
    MPI_Finalize();
    return 0;
}
