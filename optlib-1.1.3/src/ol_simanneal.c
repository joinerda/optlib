/**********************************************
*  Copyright 2007-2008, David Joiner
*
*  Simulated annealing routines
*  USAGE:   call simanneal with initial guess to minimum,
*           reference to function to be minimized,
*           and additional information
*
*           n           - number of fit parameters
*           x           - initial fit guess
*           step_factor - per-fit adjustments, 1.0 is normal
*                         step, a value <= 0.0 removes that
*                         parameter from the fit process
*                         - MAY BE SET EQUAL TO NULL -
*           T           - initial temperature
*           step        - initial step
*           step_max    - maximum step size
*           cooling     - cooling factor used to adjust change
*                         in step size and temperature each iteration
*           ITMAX       - maximum # of iterations
*           restart     - # of iterations before resetting to
*                         current best value
*           EPS         - value of step at which convergence is reached
*           func        - function to be minimized
*                         func(x,func_data) should allow
*                         for input of the parameters x and
*                         a pointer to a user provided memory
*                         structure.
*           func_data   - User provided memory structure, this
*                         allows the user if needed to pass additional
*                         information to the routine func beyond the
*                         current guess for the best fir parameters.
*                         IF UNUSED, THIS MAY BE SET EQUAL TO NULL.
*           output      - User provided routine to customize
*                         output each iteration.
*                         THIS MAY BE SET EQUAL TO NULL
*           cp_fname    - checkpoint filename, may be NULL
*
*           RECOVER     - if recovering from previous run, current
*                         iteration
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

#include "ol_simanneal.h"
#include "ol_mem.h"
#include "ol_rand_tools.h"

#define EQUILIBRIUM_TOLERANCE 0.2
#define BOLTZMAN_K 1.0

#define ANNEAL_TOGETHER 1

#ifdef USE_OPTMPI
typedef struct { double min; int loc; } MinLoc;
#endif


void simanneal(SimAnneal ** pSim,
        int n, int n_annealers,double * x,double * step_factor, double T,
        double step,double step_max, double cooling,
        int ITMAX,int restart, double EPS, int ADAPTIVE,
        double (*func)(double * x,void *),void * func_data,
        void output(double * x, void * func_data, const char * mesg),
        const char * cp_fname, int RECOVER) {
    int done=0;
    int count;
    int i;
    FILE * cp_file = NULL;

#ifdef USE_OPTMPI
    MPI_Comm_size(ebsa_mpi_comm,&ebsa_mpi_size);
    MPI_Comm_rank(ebsa_mpi_comm,&ebsa_mpi_rank);
    n_annealers = n_annealers/ebsa_mpi_size;
    if(n_annealers<1) n_annealers=1;
#endif

    if(cp_fname!=NULL) {
        cp_file = fopen(cp_fname,"w");
    }

    simanneal_init(pSim,n,n_annealers,x,step_factor,T,
        step,step_max,cooling,
        ITMAX,restart,EPS,
        func,func_data,
        output);

    if(RECOVER) (*pSim)->iter = RECOVER;

    simanneal_set_adaptive(*pSim, ADAPTIVE);
    
    count=0;
    while(!done){
        done = simanneal_iterate(*pSim,n,n_annealers,x,step_factor,step_max, cooling,
            ITMAX,restart, EPS,
            func,func_data,output);
        if(SIMANNEAL_CHECKPOINT &&
                count++%SIMANNEAL_CHECKPOINT_FREQ==0 &&
                cp_file !=NULL) {
            fprintf(cp_file,"SIMANNEAL CHECKPOINT -----------------------\n");
            fprintf(cp_file,"SIMANNEAL CHECKPOINT iter = %d\n",(*pSim)->iter);
            fprintf(cp_file,"SIMANNEAL CHECKPOINT T = %10.3e\n",(*pSim)->T);
            fprintf(cp_file,
                "SIMANNEAL CHECKPOINT step = %10.3e\n",(*pSim)->step);
            for (i=0;i<n*n_annealers;i++) {
                fprintf(cp_file,
                    "SIMANNEAL CHECKPOINT x[%d] = %10.3e\n",i,x[i]);
            }
        }
    }
    
    for(i=0;i<n;i++) {
        x[i]=((*pSim)->lowest_x)[i];
    }
    if((*pSim)->iter>=ITMAX) {
        printf("WARNING: ITMAX EXCEEDED IN SIMANNEAL\n");
        printf("         ITER = %d\n",(*pSim)->iter);
        printf("         T = %lf\n",(*pSim)->T);
        printf("         step = %lf\n",(*pSim)->step);
    }

    if(cp_file!=NULL) fclose(cp_file);
    
}

void simanneal_init(SimAnneal ** theSim,int n,int n_annealers,double * x,double * step_factor, double T,
        double step,double step_max, double cooling,
        int ITMAX,int restart, double EPS,
        double (*func)(double * x,void *),void * func_data,
        void output(double * x, void * func_data, const char * mesg)){
    int i;
    
    //initialize structure
    *theSim = (SimAnneal *)malloc(sizeof(SimAnneal));
    //allocate arrays
    (*theSim)->lowest_x = alloc_darray(n);
    (*theSim)->temp_x = alloc_darray(n);
    (*theSim)->t_step = alloc_darray(n);
    (*theSim)->new_x = alloc_darray(n);
    (*theSim)->index = alloc_iarray(n);
    (*theSim)->E = func(x,func_data);
    (*theSim)->lowest_E=(*theSim)->E;
    for(i=0;i<n;i++) (*theSim)->lowest_x[i]=x[i];
    (*theSim)->iter=0;
    (*theSim)->restart_counter=0;
    (*theSim)->sync_counter=0;
    (*theSim)->T =T;
    (*theSim)->step=step;
    (*theSim)->simanneal_step_type = SIMANNEAL_ADAPTIVE_FALSE;
}

void simanneal_finalize(SimAnneal * theSim) {
    //free arrays
    free_darray(theSim->lowest_x);
    free_darray(theSim->t_step);
    free_darray(theSim->new_x);
    free_iarray(theSim->index);
    free_darray(theSim->temp_x);
    //free structure
    free(theSim);
}

void simanneal_set_adaptive(SimAnneal * theSim, int flag) {
    if(flag) {
        theSim->simanneal_step_type = SIMANNEAL_ADAPTIVE_TRUE;
    } else {
        theSim->simanneal_step_type = SIMANNEAL_ADAPTIVE_FALSE;
    }
}

void simanneal_adjustSteps(SimAnneal * theSim,int n,int n_annealers,double * x, double * step_factor,
        double (*func)(double * x, void * func_data), void * func_data){
    
    //OK, start with each guess, take a positive step, a negative step, and
    // compute the difference in energy--also graceful errors for negative,zero setp_factor, NULL stepfactor needed
    double * test_x;
    double test_E;
    double * delta_E;
    int i,j,i_min;
    double delta_Emin;

    
    test_x = alloc_darray(n);
    delta_E = alloc_darray(n);

    delta_Emin=-1.0;
    i_min=0;
    for(i=0; i<n;i++) {
        if(step_factor[i]>=0.0){
            step_factor[i]=1.0;
            for(j=0;j<n;j++) {
                test_x[j]=theSim->lowest_x[j];
            }
            test_x[i]=theSim->lowest_x[i]+step_factor[i]*theSim->step;
            test_E = func(test_x,func_data);
            test_x[i]=theSim->lowest_x[i]-step_factor[i]*theSim->step;
            delta_E[i]=fabs(func(test_x,func_data)-test_E);
            if((delta_E[i]<delta_Emin||delta_Emin<0.0)&&delta_E[i]>0.0) {
                i_min=i;
                delta_Emin=delta_E[i];
            }
        } else {
            delta_E[i]=-1.0;
        }
    }
    
    for(i=0;i<n;i++) {
        if(delta_E[i]>0.0&&i!=i_min){
            step_factor[i]*=max(delta_Emin/delta_E[i],0.1);
        }
    }
    
    free_darray(test_x);
    free_darray(delta_E);
}

int simanneal_iterate(SimAnneal * theSim,int n,int n_annealers,double * x,double * step_factor,
        double step_max, double cooling,
        int ITMAX,int restart, double EPS,
        double (*func)(double * x,void *),void * func_data,
        void output(double * x, void * func_data, const char * mesg)) {
    int i,j,k;
    int step_succeeded;
    int ITMIN=3;
    int done;
    int i_inner,inner_max,trials_per_temp;
    int attempts, good_attempts, bad_attempts;
    int do_restart;
    int iter_change,counter_change;
#ifdef USE_OPTMPI
    int attempts_buf,good_attempts_buf,bad_attempts_buf,iter_buf,counter_buf;
    MinLoc minloc,minloc_buf;
    int done_buf;
#endif

    done=0;
    trials_per_temp=200;
    attempts=0;
    good_attempts=0;
    bad_attempts=0;
    do_restart=0;
#ifdef USE_OPTMPI
    attempts_buf=0;
    good_attempts_buf=0;
    bad_attempts_buf=0;
#endif
    iter_change=0;
    counter_change=0;
#ifdef USE_OPTMPI
    inner_max=(int)((double)trials_per_temp/(double)n_annealers/
        (double)ebsa_mpi_size);
    if(inner_max<1) {
        printf("WARNING: INNER_MAX less than 1 in simanneal\n");
        printf("        %d   %d   %d\n",trials_per_temp, n_annealers, ebsa_mpi_size);
        inner_max=1;
    }
#else
    inner_max=(int)((double)trials_per_temp/(double)n_annealers);
#endif
    // split here across processes, can this be done SIMD?
    //   every process does a portion of the annealers and shares
    //   number of sucesses? At end of turn, each shares the best
    //   current answer, at each restart, everyone is set to the
    //   global best answer.
    //   At each checkpoint, the entire solution is written to
    //   the checkpoint file.
    //   Each process initializes its own annealers.
    for (j=0;j<n_annealers;j++) {
        for(i_inner=0;i_inner<inner_max;i_inner++) {
            for(k=0;k<n;k++) {
                theSim->temp_x[k]=x[n*j+k];
            }
            step_succeeded = simanneal_step(theSim,n,(theSim->temp_x),step_factor,
                &(theSim->E),theSim->T,theSim->step,func,func_data,output);
            if(step_succeeded)
                for(k=0;k<n;k++) x[n*j+k]=theSim->temp_x[k];
            attempts++;
            if(step_succeeded) {
                good_attempts++;
            } else {
                bad_attempts++;
            }
            iter_change++;
            counter_change++;
            if(theSim->E<theSim->lowest_E) {
                theSim->lowest_E=theSim->E;
                for(i=0;i<n;i++) theSim->lowest_x[i]=x[j*n+i];
            }
        }
    }
#ifdef USE_OPTMPI
    MPI_Allreduce(&attempts,&attempts_buf,1,MPI_INT,MPI_SUM,ebsa_mpi_comm);
    MPI_Allreduce(&good_attempts,&good_attempts_buf,1,
        MPI_INT,MPI_SUM,ebsa_mpi_comm);
    MPI_Allreduce(&bad_attempts,&bad_attempts_buf,1,
        MPI_INT,MPI_SUM,ebsa_mpi_comm);
    MPI_Allreduce(&iter_change,&iter_buf,1,
        MPI_INT,MPI_SUM,ebsa_mpi_comm);
    MPI_Allreduce(&counter_change,&counter_buf,1,
        MPI_INT,MPI_SUM,ebsa_mpi_comm);
    attempts = attempts_buf;
    good_attempts = good_attempts_buf;
    bad_attempts = bad_attempts_buf;
    theSim->iter+=iter_buf;
    theSim->restart_counter+=counter_buf;
#else 
    theSim->iter+=iter_change;
    theSim->restart_counter+=counter_change;
#endif
    if(theSim->restart_counter>restart) {
        do_restart=1;
        theSim->restart_counter=0;
    }
    if(do_restart&&restart>0) {
        for(j=0;j<n_annealers;j++) {
            for(i=0;i<n;i++) {
                theSim->E=theSim->lowest_E;
                for(i=0;i<n;i++) x[j*n+i]=theSim->lowest_x[i];
                theSim->T *= 2.0;
                theSim->step *= 10.0;
            }
        }
    }
#ifdef USE_OPTMPI
    MPI_Bcast(&(theSim->sync_counter),1,MPI_INT,0,ebsa_mpi_comm);
#endif
    if(ANNEAL_TOGETHER&&theSim->sync_counter%100==0) {
        int jkeep=0;
        int jmin=0;
        double Emin = func(&(x[0*n]),func_data);
        double Etest;
        double Ekeep;
        double dtest;
        double tbuffer;
        jkeep = irand(0,n_annealers-1);
        for(j=1;j<n_annealers;j++) {
            Etest = func(&(x[j*n]),func_data);
            if(j==jkeep) Ekeep=Etest;
            if(Etest<Emin) {
                jmin=j;
                Emin=Etest;
            }
        }
        if(drand(0,1)>exp((Etest-Ekeep)/(BOLTZMAN_K*theSim->T))) {
            jkeep=jmin;
            Ekeep=Emin;
        } else {
            Ekeep=Etest;
        }
        for(j=0;j<n_annealers;j++) {
            if(j!=jkeep) {
                for(i=0;i<n;i++) {
                    x[j*n+i]=x[jkeep*n+i];
                }
            }
        }
#ifdef USE_OPTMPI
        minloc.min=Ekeep;
        minloc.loc=ebsa_mpi_rank;
        MPI_Allreduce(&minloc,&minloc_buf,1,
            MPI_DOUBLE_INT,MPI_MINLOC,ebsa_mpi_comm);
        minloc.min = minloc_buf.min;
        minloc.loc = minloc_buf.loc;
        if(ebsa_mpi_rank==0) jkeep = irand(0,ebsa_mpi_size);
        MPI_Bcast(&jkeep,1,MPI_INT,0,ebsa_mpi_comm);
        MPI_Bcast(&Ekeep,1,MPI_DOUBLE,jkeep,ebsa_mpi_comm);
        if(ebsa_mpi_rank==0) dtest = drand(0,1);
        MPI_Bcast(&dtest,1,MPI_DOUBLE,0,ebsa_mpi_comm);
        MPI_Allreduce(&(theSim->T),&tbuffer,1,MPI_DOUBLE,MPI_MIN,ebsa_mpi_comm);
        if(dtest>exp((minloc.min-Ekeep)/(BOLTZMAN_K*tbuffer))) {
            jkeep=minloc_buf.loc;
        } 
        MPI_Bcast(x,n,MPI_DOUBLE,jkeep,ebsa_mpi_comm);
        for(j=1;j<n_annealers;j++) {
             for(i=0;i<n;i++) {
                 x[j*n+i]=x[i];
             }
        }
#endif
    }
    theSim->sync_counter++;
    if(SIMANNEAL_VERBOSE||SIMANNEAL_CHECKPOINT) {
#ifdef USE_OPTMPI
        if((ebsa_mpi_rank==0&&ANNEAL_TOGETHER)||(!ANNEAL_TOGETHER)) {
#endif
            printf("SIMANNEAL: ITERATION %d T %10.3e STEP %10.3e \n",
                theSim->iter,theSim->T,theSim->step);
            if (output != NULL) {
                output(x,func_data,       "SIMANNEAL: CURRENT X ");
                output(theSim->lowest_x,func_data,"SIMANNEAL: LOWEST  X ");
            } else {
                print_darray(n*n_annealers,x,       "SIMANNEAL: CURRENT X ");
                print_darray(n,theSim->lowest_x,"SIMANNEAL: LOWEST  X ");
            }
            printf("SIMANNEAL: ENERGY, LOWEST ENERGY  %13.6e %13.6e\n",
                theSim->E,theSim->lowest_E);
#ifdef USE_OPTMPI
        }
#endif
    }
    if(fabs((double)bad_attempts/(double)attempts-0.5)<EQUILIBRIUM_TOLERANCE) {
        theSim->T *= cooling;
        if(theSim->simanneal_step_type==SIMANNEAL_ADAPTIVE_TRUE&&
                step_factor!=NULL)
            simanneal_adjustSteps(theSim,n,n_annealers,x,step_factor,func,func_data);
    } else {
        if(bad_attempts>good_attempts) { 
            theSim->step *= (1.0-0.1/(double)n_annealers);
        } else {
            theSim->step *= (1.0+0.1/10.0/(double)n_annealers);
        }
    }
    if(theSim->step>step_max) {
        theSim->step=step_max/2.0;
        theSim->T*=cooling;
    }
    
    // determine convergence
    done = 1;
    if(step_factor==NULL) {
        // no step factor used, try relative first, then absolute.
        for (i=0;i<n*n_annealers;i++) {
            if (fabs(x[i])>EPS) {
                if (fabs(theSim->step/x[i])>EPS) {
                    done=0;
                }
            } else {
                if (theSim->step>EPS) {
                    done=0;
                }
            }
        }
    } else {
        // step factor used,try relative first,then absolute,for absolute,use
        // larger of step size and scaled step size to determine convergence
        for (i=0;i<n*n_annealers;i++) {
            if (step_factor[i%n]>0.0) {
                if (fabs(x[i])>EPS) {
                    if (fabs(theSim->step/x[i])>EPS) {
                        done=0;
                    }
                } else {
                    if (theSim->simanneal_step_type==SIMANNEAL_ADAPTIVE_TRUE) {
                        if(max(theSim->step,
                               fabs(step_factor[i%n_annealers]*theSim->step))>EPS) {
                            done=0;
                        }
                    } else {
                        if(theSim->step>EPS) {
                            done=0;
                        }
                    }
                }
            }
        }
    }
    if(theSim->iter<ITMIN) {
        done=0;
    }
    if(theSim->iter>ITMAX) {
        done=1;
    }
#ifdef USE_OPTMPI
    MPI_Allreduce(&done,&done_buf,1,
        MPI_INT,MPI_PROD,ebsa_mpi_comm);
    done = done_buf;
#endif

    return done;
}

int simanneal_step(SimAnneal * theSim, int n_inp, double *x,double * step_factor,
        double *E, double T,double delta,
        double (*func)(double *x,void *),void * func_data,
        void output(double * x, void * func_data, const char * mesg)) {
    double new_E;
    int i,m;
    int accepted;
    int step_method;
    int step_chosen;
    int n;
    double Etest;

    double change;


    for(i=0;i<n_inp;i++) theSim->new_x[i]=0.0;
    if(step_factor!=NULL) {
        n=0;
        for(i=0;i<n_inp;i++) {
            if(step_factor[i]>0.0) {
                theSim->index[n]=i;
                n++;
            }
        }
        if(n==0) {
            printf("WARNING: All step factors <= 0.0 in simanneal_step\n");
            exit(0);
        }
    } else {
        n = n_inp;
        for(i=0;i<n;i++) theSim->index[i]=i;
    }

    accepted=0;
    // choose type of change: pick a variable to vary or move in "all space"
    step_chosen = 0;
    while(!step_chosen) {
        step_method = rand()%(3*n);
        if(step_method<n&&SIMANNEAL_STEP_METHOD_1) {
            if(step_factor==NULL) {
                step_chosen=1;
            } else if(step_factor[theSim->index[step_method]]>0.0) {
                step_chosen=1;
            }
            if(step_chosen) {
                for (i=0;i<n_inp;i++) theSim->new_x[i]=x[i];
                change = drand(-1,1)*delta;
                if(step_factor!=NULL) change *= step_factor[theSim->index[step_method]];
                theSim->new_x[theSim->index[step_method]]+=change;
            }
        } else {
            if(step_method<2*n&&SIMANNEAL_STEP_METHOD_1) {
                m = irand(1,max(n/2,1));
            } else {
                m = irand(1,n);
            }
            random_direction_subset(n,m,theSim->t_step);
            for (i=0;i<n_inp;i++) theSim->new_x[i]=x[i];
            for (i=0;i<n;i++) {
                if(step_factor!=NULL) {
                    theSim->new_x[theSim->index[i]]+=theSim->t_step[i]*delta*step_factor[theSim->index[i]];
                } else {
                    theSim->new_x[theSim->index[i]]+=theSim->t_step[i]*delta;
                }
            }
            step_chosen=1;
        }
    }
    if(SIMANNEAL_VERBOSE) {
/*
        if (output != NULL) {
            output(theSim->new_x,func_data,  "SIMANNEAL: ATTEMPT X ");
        } else {
            print_darray(n_inp,theSim->new_x,"SIMANNEAL: ATTEMPT X ");
        }
*/
    }
    new_E = func(theSim->new_x,func_data); 
    Etest=*E;
    if((drand(0,1))<(exp((Etest-new_E)/(BOLTZMAN_K*T)))) {
        accepted=1;
        for(i=0;i<n_inp;i++) x[i]=theSim->new_x[i];
        *E=new_E;
    }
    

    return accepted;
}
