#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <optlib.h>

// determination of vertical miss distance squared,
//  assuming target is 100 m away and 50 m up.
double target(double *guess, void * func_data) {
    double v = guess[0];
    double theta = guess[1];

    double x,y,vx,vy,t,dt;
    x = 0.0;
    y = 0.0;
    vx = v*cos(theta);
    vy = v*sin(theta);
    dt = 100.0/v/100.0;
    if(dt>0.01) dt=0.01;
    t = 0;
    if(theta<0.0||theta>M_PI/2) {
        return 10000.0+theta*theta;
    }
    if(vx<0.0) { // negative vx penalty
        return 10000.0+vx*vx;
    }
 //   if(v>40.0) { // large v penalty
 //       return 10000.0+v*v;
 //   }
    while(x<100.0) {
        x += vx*dt;
        y += vy*dt;
        vy += -9.8*dt;
        t += dt;
        if(y<0.0) {
            return 10000.0+y*y; //landing too early penalty
        }
    }
    return (y-50.0)*(y-50.0) ; //squared to ensure minimum and
                              // not saddle point;
}

int main(int argc, char ** argv) {
    double guess[2];
    int i=0;
    int seed_offset=0;
    OptLibOpts * theOpts;

    // always seed your stochastic models
    if(argc>1) sscanf(argv[1],"%d",&seed_offset);
    seed_by_time(seed_offset);
    
    theOpts = OPTLIB_CreateOpts();
    printf("v theta e\n");
    for (i=0;i<20;i++) {
        // inintialize guess
        guess[0]=drand(0,100.0);
        guess[1]=drand(0,M_PI/2.0);
    
        // run OPTLIB_Minimize with defaults
        theOpts->ga_itmax=200;
        theOpts->ga_itmin=10;
        theOpts->ga_dominancefactor=0.4;
        OPTLIB_Minimize(2,guess,&target,NULL,theOpts);
    
        // output
        printf("%lf %lf %lf\n",
            guess[0],guess[1],
            target(guess,NULL));
    }
    OPTLIB_DestroyOpts(theOpts);

   
    return 0;
}
