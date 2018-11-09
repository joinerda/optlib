#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <optlib.h>

#define MAX_LINE_LENGTH 240

// data structure used to pass data to chi_squared routine
typedef struct {
    int n;
    double * x;
    double * y;
} data;

double chi_squared(double *x, void * func_data) {
    double a=x[0];
    double b=x[1];
    data * theData = (data *)func_data;
    double sum;
    int i;

    // loop over data and calculate chi-squared assuming model
    //  of a*exp(-b*x*x)
    sum=0.0;
    for(i=0;i<theData->n;i++) {
        double f = a*exp(-b*theData->x[i]*theData->x[i]);
        sum += (f-theData->y[i])*(f-theData->y[i])/fabs(theData->y[i]);
    }

    return sum;
    

}

int main(int argc, char ** argv) {
    FILE * infile;
    char infile_name[MAX_LINE_LENGTH];
    char line[MAX_LINE_LENGTH];
    data theData;
    double guess[2];
    int i=0;


    // always seed your stochastic models
    seed_by_time(0);
    
    // read in data file
    // file format should be number of data items on first line,
    //   x,y, values on successive lines
    sprintf(infile_name,"data.txt");
    infile = fopen(infile_name,"r");
    theData.x=NULL;
    theData.y=NULL;
    theData.n=0;
    while(fgets(line,MAX_LINE_LENGTH,infile)!=NULL) {
        if(theData.n==0) {
            sscanf(line,"%d",&(theData.n));
            if(theData.n<1) {
                printf("ERROR: value of n=%d not allowable\n",theData.n);
                exit(0);
            }
            theData.x = (double *)malloc(sizeof(double)*theData.n);
            theData.y = (double *)malloc(sizeof(double)*theData.n);
        } else {
            sscanf(line,"%lf %lf",&(theData.x[i]),&(theData.y[i]));
            i++;
        }
    }
    if(i!=theData.n) {
        printf("ERROR: value of i=%d not equal to n=%d\n",i,theData.n);
        exit(0);
    }
    fclose(infile);

    // inintialize guess
    guess[0]=1.0;
    guess[1]=1.0;

    // run OPTLIB_Minimize with defaults
    OPTLIB_Minimize(2,guess,&chi_squared,&theData,NULL);

    // output
    printf("SOLUTION a[%lf] b[%lf] chi-squared[%lf]\n",guess[0],guess[1],
        chi_squared(guess,&theData));

    // free data
    if(theData.x!=NULL) free(theData.x);
    if(theData.y!=NULL) free(theData.y);
   
    return 0;
}
