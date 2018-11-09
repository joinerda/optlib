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

#include <stdio.h>
#include <math.h>
#include "ol_mem.h"
#include "ol_rand_tools.h"
#include "ol_genetic_alg.h"

#define TINY 1.0e-30

#define GA_VERBOSE_ALG_RESULTS 0
#define GA_VERBOSE_POOLS_RESULTS 0

void genetic_mutation(int n, double *mutation, double *original,
    int n_step, double * step_factor) {
    int j;
    double range;

    for(j = 0;j<n;j++) {
        if(step_factor==NULL) {
            range = fabs(original[j])*2.0;
            mutation[j] = original[j] + drand(-range,range);
        } else {
            if(n_step==1) {
                range = 2.0*step_factor[0]*fabs(original[j]);
                if(range<TINY) range=step_factor[0];
                mutation[j] = original[j] + drand(-range,range);
            } else if(n_step==n) {
                range = 2.0*max(step_factor[j],0.0)*fabs(original[j]);
                if(range<TINY) range=step_factor[j];
                mutation[j] = original[j] + drand(-range,range);
            } else {
                printf("ERROR: In genetic_mutation, n_step must ");
                printf("equal 1 or n\n");
            }
        }
    }
}


void quicksort(int n, double * x, int * index) {
    double pivot;
    int i;
    int * less;
    int i_pivot,i_less,i_greater;
    int * greater;

    if(n>1) {
        less = alloc_iarray(n);
        greater = alloc_iarray(n);
        i_less=0;
        i_greater=0;
        i_pivot = index[n/2];
        pivot = x[i_pivot];
        for (i=0;i<n;i++) {
            if(index[i]!=i_pivot) {
                if(x[index[i]]<pivot) {
                    less[i_less++]=index[i];
                } else {
                    greater[i_greater++]=index[i];
                }
            }
        }
        for(i=0;i<i_less;i++) {
            index[i]=less[i];
        }
        index[i_less]=i_pivot;
        for(i=0;i<i_greater;i++) {
            index[n-1-i]=greater[i];
        }
        quicksort(i_less,x,index);
        quicksort(i_greater,x,&(index[n-i_greater]));
        free_iarray(less);
        free_iarray(greater);
    }

    return;

}

void genetic_combination(int n,double * child, double * mother,
        double * father,
        int n_step,double * step_factor,double mutation_factor,
        double dominance_factor,double * range) {
    double test,average,sigma,sigmamin;
    double chosen;
    int j;

    for(j=0;j<n;j++) {
        average = 0.5*(mother[j]+father[j]);
        test=drand(0,1);
        if(test<mutation_factor) {
            sigma = max(average,10.0*range[j]);
            chosen = average;
        } else if(test>1.0-dominance_factor) {
            if(drand(0,1)<0.5) chosen = mother[j];
            else chosen = father[j];
            sigma = min(chosen,0.1*range[j]);
        } else {
            sigma = range[j];
            chosen = average;
        }
        if(step_factor!=NULL) {
            if(n_step==1) {
                sigma *=step_factor[0];
            } else if (n_step==n) {
                sigma *= max(step_factor[j],0.0);
            } else {
                printf("ERROR: In genetic_combination, ");
                printf("n_step must = 1 or n\n");
            }
        }
        child[j] = drand_norm(chosen,sigma,0.01);
    }
}

void genetic_mix(int n,int n_pop, int n_keep,int itmax,int itmin,
        double epsilon,
        int n_step, double * step_factor, double mutation_factor,
        double dominance_factor, double ** population,
        void output(double * x, void * func_data, const char * mesg),
        fw_struct * func_struct) {
    double * energy;
    double * average;
    double * average_squared;
    double * deviation;
    double max_deviation;
    double ** temp_population;
    double * temp_energy;
    double * range;
    double * range_max;
    double * range_min;
    int * index;
    int i,j;
    int mother,father;
    int done;
    int count;
    extern int printbest;

    //seed_by_time(0);

    energy = alloc_darray(n_pop);
    average = alloc_darray(n);
    average_squared = alloc_darray(n);
    deviation = alloc_darray(n);
    temp_energy = alloc_darray(n_keep);
    range = alloc_darray(n);
    range_max = alloc_darray(n);
    range_min = alloc_darray(n);
    temp_population = alloc_dmatrix(n_pop,n);
    index = alloc_iarray(n_pop);

    done = 0;

    count=0;
    while(!done) {
        for(i=0;i<n_pop;i++) {
            //energy[i] = func_struct->func(population[i],func_data);
            if(GA_MPI_VERBOSE&&i<5) {
                if(output!=NULL) {
                    output(population[i],
                        func_struct->func_data,"DEBUG: In GENETIC_MIX\t");
                } else {
                    print_darray(n,population[i],"DEBUG: In GENETIC_MIX\t");
                }
            }
        }
        fw_simd_func(population,energy,n_pop,func_struct);
        for (i=0;i<n_pop;i++) index[i]=i;
        quicksort(n_pop,energy,index);
        
        for (j=0;j<n;j++) {
            range_max[j] = population[index[0]][j];
            range_min[j] = population[index[0]][j];
        }
        for (i=0;i<n_keep;i++) {
            temp_energy[i]=energy[index[i]];
            for(j=0;j<n;j++) {
                temp_population[i][j]=population[index[i]][j];
                if(population[index[i]][j]>range_max[j])
                    range_max[j]=population[index[i]][j];
                if(population[index[i]][j]<range_min[j])
                    range_min[j]=population[index[i]][j];
            }
        }
        for (j=0;j<n;j++) {
            range[j] = max(epsilon,range_max[j]-range_min[j]);
        }
        for (i=0;i<n_keep;i++) {
            energy[i]=temp_energy[i];
            for(j=0;j<n;j++) {
                population[i][j]=temp_population[i][j];
            }
        }
        for (i=n_keep;i<n_pop;i++) {
            mother = irand(0,n_keep);
            do {
                father = irand(0,n_keep);
            } while (father!=mother);
            genetic_combination(n,population[i],population[mother],
                    population[father],n_step,step_factor,mutation_factor,
                    dominance_factor,range);
        }

        count++;

        max_deviation=0.0;
        if(GA_MPI_VERBOSE) printf("DEBUG count = %d \n",count);
        j=0;
        average[j] = 0.0;
        average_squared[j] = 0.0;
        for(i=0;i<n_keep;i++) {
            average[j] += energy[i];
            average_squared[j] += pow(energy[i],2.0);
        }
        average[j] /= (double)n_keep;
        average_squared[j] /= (double)n_keep;
        deviation[j] = pow(fabs(-average[j]*average[j]+
            average_squared[j]),0.5)*min(1.0,1.0/fabs(average[j]));
        if(GA_MPI_VERBOSE) printf("DEBUG deviation = %d %g\n",j,deviation[j]);
        if(deviation[j]>max_deviation) {
            max_deviation=deviation[j];
        }
        /*
        for(j=0;j<n;j++) {
            average[j] = 0.0;
            average_squared[j] = 0.0;
            for(i=0;i<n_keep;i++) {
                average[j] += population[i][j];
                average_squared[j] += pow(population[i][j],2.0);
            }
            average[j] /= (double)n_keep;
            average_squared[j] /= (double)n_keep;
            deviation[j] = pow(fabs(-average[j]*average[j]+
                average_squared[j]),0.5)*min(1.0,1.0/fabs(average[j]));
            if(GA_MPI_VERBOSE)
                printf("DEBUG deviation = %d %g\n",j,deviation[j]);
            if(deviation[j]>max_deviation) {
                max_deviation=deviation[j];
            }
        }
        */
        if(max_deviation<epsilon&&count>itmin) done=1;
        if(count>itmax) done=2;
    }

    if(done==2 && epsilon>0.0) {
        printf("WARNING: itmax exceeded in genetic_mix()\n");
    }


    free_dmatrix(temp_population);
    free_darray(energy);
    free_darray(average);
    free_darray(average_squared);
    free_darray(deviation);
    free_darray(temp_energy);
    free_darray(range);
    free_darray(range_max);
    free_darray(range_min);
    free_iarray(index);
}

void genetic_alg_serial(int n, int n_pop, int n_keep, int itmax,
    int itmin, double epsilon, int n_step, double * step_factor,
    double mutation_factor, double dominance_factor, double * guess,
    void output(double * x, void * func_data, const char * mesg),
    double (*func)(double *, void *), void * func_data) {

    fw_struct * func_struct = (fw_struct *)malloc(sizeof(func_struct));
    func_struct->func = func;
    func_struct->func_data = func_data;
    func_struct->n = n;
    func_struct->size = 1;
    func_struct->rank = 0;
    
    genetic_alg(n,n_pop,n_keep,itmax,itmin,epsilon,
        n_step,step_factor,mutation_factor,
        dominance_factor, guess, output, func_struct);

    free(func_struct);
}

void genetic_alg(int n, int n_pop, int n_keep,int itmax, int itmin, double epsilon,
        int n_step, double *step_factor, double mutation_factor,
        double dominance_factor,
        double * guess,
        void output(double * x, void * func_data, const char * mesg),
        fw_struct * func_struct) {
    double ** population;
    int i,j;

    population = alloc_dmatrix(n_pop,n);

    for (j=0;j<n;j++) {
        population[0][j] = guess[j];
    }
    for (i=1;i<n_pop;i++) {
        genetic_mutation(n,population[i],guess,n_step,step_factor);
    }
    genetic_mix(n,n_pop,n_keep,itmax,itmin,epsilon,n_step,step_factor,
        mutation_factor,dominance_factor,
        population,output,func_struct);
    for(j=0;j<n;j++) guess[j]=population[0][j];

    if(GA_VERBOSE_ALG_RESULTS) {
        for(i=0;i<n_pop;i++) {
            print_darrayt(n,population[i],"GA_RESULT ");
        }
    }

    free_dmatrix(population);
}



void genetic_pools_serial(int n, int n_pools, int n_pop, int n_keep,
        int itmax, int itmin,
        double epsilon, int n_step, double * step_factor,
        double mutation_factor,
        double dominance_factor,
        double * guess, 
        void output(double * x, void * func_data, const char * mesg),
        double (*func)(double *, void *), void * func_data) {

    fw_struct * func_struct = (fw_struct *)malloc(sizeof(fw_struct));
    func_struct->func = func;
    func_struct->func_data = func_data;
    func_struct->size = 1;
    func_struct->rank = 0;
    func_struct->n = n;

    genetic_pools(n,n_pools,n_pop,n_keep, itmax,itmin,
        epsilon,n_step,step_factor, mutation_factor,
        dominance_factor, guess, output, func_struct);

    free(func_struct);
}

void genetic_pools(int n, int n_pools, int n_pop, int n_keep,
        int itmax, int itmin,
        double epsilon, int n_step, double * step_factor,
        double mutation_factor,
        double dominance_factor,
        double * guess, 
        void output(double * x, void * func_data, const char * mesg),
        fw_struct * func_struct) {
    double ** population;
    int i,j;

    population = alloc_dmatrix(n_pop,n);

    for (i=0;i<n_pools;i++) {
        if(i==0) {
            for(j=0;j<n;j++) {
                population[i][j]=guess[j];
            }
        } else {
            genetic_mutation(n,population[i],guess,n_step,step_factor);
        }
        genetic_alg(n,n_pop,n_keep,itmax,itmin,epsilon,n_step,step_factor,
            mutation_factor, dominance_factor,
            population[i],output,func_struct);
        if(GA_MPI_VERBOSE) print_darray(n,population[i],"DEBUG -- \t");
    }
    for (i=n_pools;i<n_pop;i++) {
        for(j = 0;j<n;j++) {
            population[i][j]=population[i%n_pools][j];
        }
        if(GA_MPI_VERBOSE)print_darray(n,population[i],"DEBUG -- \t");
    }
    genetic_mix(n,n_pop,n_keep,itmax,itmin,epsilon,n_step,step_factor,
        mutation_factor,dominance_factor,
        population,output,func_struct);
    for (j=0;j<n;j++) guess[j]=population[0][j];

    if(GA_VERBOSE_POOLS_RESULTS) {
        for(i=0;i<n_pop;i++) {
            print_darray(n,population[i],"GA_POOLS_RESULT ");
        }
    }

    free_dmatrix(population);
}

