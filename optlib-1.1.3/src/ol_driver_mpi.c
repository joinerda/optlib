#include <optlib_mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include "func.h"
#ifdef USE_OPTMPI
#include <mpi.h>
#endif

int main(int argc, char ** argv) {
    OptLibOpts * theOpts;
    int n=2;
    double  *guess;
    int i;
    int size=1;
    int rank=0;
    func_data_struct func_data;
    func_data.fcount=0;
    int buffer;

#ifdef USE_OPTMPI
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&size);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    printf("MPI Driver version, message from rank %d/%d\n",rank,size);
#else
    printf("WARNING: MPI Driver version, ");
    printf("but USE_OPTMPI not set at compile time\n");
#endif
    if(rank==0) {
        printf("OPTLIB DRIVER version %s.%s\n",optlib_VERSION_MAJOR,optlib_VERSION_MINOR);
    }

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
    for (i=0;i<3;i++) {
        func_data.fcount=0;
        func_data.method=i;
        theOpts->method = i;

        guess[0]=drand(0,1);
        guess[1]=drand(0,1);
        OPTLIB_Minimize(n,guess,&func,(void *)&func_data,theOpts);

#ifdef USE_OPTMPI
        MPI_Reduce(&(func_data.fcount),&buffer,1,MPI_INT,
            MPI_SUM,0,MPI_COMM_WORLD);
        func_data.fcount=buffer;
#endif
        if(rank==0) {
            if (theOpts->method==OPTLIB_METHOD_GA) printf("METHOD GA\n");
            else if (theOpts->method==OPTLIB_METHOD_SA) printf("METHOD SA\n");
            else if (theOpts->method==OPTLIB_METHOD_MCP) printf("METHOD MCP\n");
            else if (theOpts->method==OPTLIB_METHOD_POWELL)
                printf("METHOD POWELL\n");
            printf("SOLUTION %g %g\n",guess[0], guess[1]);
            printf("TOTAL FUNCTION CALLS %d\n",func_data.fcount);
        }
    }

    OPTLIB_DestroyOpts(theOpts);
    free(guess);
#ifdef USE_OPTMPI
    MPI_Finalize();
#endif
    return 0;
}
