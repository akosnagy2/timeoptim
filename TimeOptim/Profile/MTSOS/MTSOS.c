#ifdef MATLAB_MEX_FILE
    #include "mex.h"
    #define malloc mxMalloc
    #define free mxFree
    #define realloc mxRealloc
    #define calloc mxCalloc
    #include "blas.h" //to ensure that the blas library is included append the argument -lmwblas to the compile
#endif
#ifndef MATLAB_MEX_FILE
	#include <cblas.h>
#endif
#include "MTSOS.h"
#include <string.h>
#include <math.h>
#include <time.h>
        
#ifndef INFINITY
    #define INFINITY HUGE_VAL
#endif

#ifdef __cplusplus
}
#endif

int sgn(double a)
{
	if (a >= 0)
		return 1;
	else
		return -1;
}

void writeOutPoint(FILE *fd, double* data);

int so_MTSOS(problem_params* p_params, algorithm_flags* flags, algorithm_params* a_params, optional_params* o_params, double** p_b, double** p_u, double** p_v, double* p_fx, double** p_timers, int* p_iterations){
//the primary function that performs optimization on the provided vector S
    double kappa, alpha, beta, mu, epsilon, sum, residual, residual_new, Lambda, time, b_0;
    double* S_middle, *b, *a, *x, *u, *v, *S;
	double *S_prime, *S_dprime, *S_e_prime, *S_e_dprime, *S_e_prime_s, dtheta, *resid, *data, *Abig_v, *bbig, *dx, *G, *indexH_v, *indexA21_v, *b_Ax, *timers, *variables;
    int* indexH_i, *indexH_j, *indexA21_i, *indexA21_j;
    int iterations, solved, feasible, display;
    int S_length, State_size, index_length, U_size, Hindex_length, A21index_length, variables_length;
    int i, j, k, status, totime, orderOkay;
    int* Abig_i, *Abig_j;
    int MAX_ITERATIONS, flags_kappa;
    double RESIDUAL_THRESHHOLD, fx;
    cs *csA21, Abig, *Abig_c, *Hc;
    css* lu_S;
    clock_t start, start2, diff;
    
    status = 1;//flag that says everything is alright

    if(flags->timer == 1){
        start = clock();
        totime = 1;
        if(*p_timers == NULL)
            timers = calloc(19, sizeof(double));
        else
            timers = realloc(*p_timers, 19*sizeof(double));
        if(timers==NULL){
            #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
            #endif
            #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
            #endif
                    
            status = 0;
            totime = 0;
        }
        *p_timers = timers;
    }else{
        totime = 0;
        timers = NULL;
    }
    
    //determine kappa, for the barrier function;
    if(a_params->kappa <= 0)
        kappa = .1;
    else
        kappa = a_params->kappa;
    
    //determine the backtracking line search parameter alpha if set
    if(a_params->alpha <= 0 || a_params->alpha >= .5)
        alpha = .01;
    else
        alpha = a_params->alpha;
    
    //determine the backtracking line searh parameter beta if set
    if(a_params->beta <= 0 || a_params->beta >= 1)
        beta = .75;
    else
        beta = a_params->beta;

    //determine the espilon, accuracy term
    if(a_params->epsilon <= 0)
        epsilon = .1;
    else
        epsilon = a_params->epsilon;

    //set the maximum number of iterations
    if(a_params->MAX_ITER <= 0)
        MAX_ITERATIONS = 100;
    else
        MAX_ITERATIONS = a_params->MAX_ITER;
    
    //see if any optional variables have been set
    if(o_params == NULL || o_params->variables == NULL){
        variables = NULL;
        variables_length = 0;
    }else{
        variables = o_params->variables;
        variables_length = o_params->variables_length;
    }
    
    mu = 10;//the increase term for the barrier funciton
    solved = 0;//has a solution been found yet;
    feasible = 0;//is the current solution feasible
    
    fx = -1;//initialize fx to -1 to show that no feasible solution was found
    Abig_v = NULL;//initialize Abig to NULL to show that space has not been allocated yet;
    bbig = NULL;//iniitalize bbig to NULL to show the same has not been allocated
    iterations = 0;//you have not yet performed any iterations

    RESIDUAL_THRESHHOLD = 1e-6;//error threshhold required for the problem to be feasible.
    orderOkay = 0;//flag for detecting anamolies in factorization requiring refactorization
    lu_S = NULL;//for storing factorization
    residual = INFINITY;//the residual has not yet been calculated
    
    U_size = p_params->U_size;
    S_length = p_params->S_length;
    if(S_length == 1){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("Insufficient trajectory, 1 point, provided");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "Insufficient trajectory, 1 point, provided\n");
        #endif
        status = 0;
    }
    State_size = p_params->State_size;
    
    //defaults to dynamic kappa (flags_kappa = 1);
    if(flags->kappa != 0)
        flags_kappa = 1;
    else
        flags_kappa = 0;

    //defaults to no display (display = 1)
    if(flags->display != 1)
        display = 0;
    else
        display = 1;
        
    if(display && status){
        #ifdef MATLAB_MEX_FILE
            mexPrintf("line search parameter alpha is %f\nline search parameter beta is %f\nbarrier parameter kappa is %f\naccuracy parameter epsilon is %f\n", alpha, beta, kappa, epsilon);
            mexPrintf("maximum number of iterations is %d\n", MAX_ITERATIONS);
        #endif
        #ifndef MATLAB_MEX_FILE
            printf("line search parameter alpha is %f\nline search parameter beta is %f\nbarrier parameter kappa is %f\naccuracy parameter epsilon is %f\n", alpha, beta, kappa, epsilon);
            printf("maximum number of iterations is %d\n", MAX_ITERATIONS);
        #endif
    }
        
    dtheta = 1.0/(S_length-1);//calcuate the length of the interval
    //calculate S_middle;
    S = p_params->S;
    S_middle = calloc(State_size*(S_length-1),sizeof(double));
    S_prime = calloc(State_size*S_length-1,sizeof(double));
    S_dprime = calloc((S_length-1)*State_size, sizeof(double));
	S_e_prime = calloc(S_length*State_size, sizeof(double));
	S_e_dprime = calloc(S_length*State_size, sizeof(double));
	S_e_prime_s = calloc(S_length*State_size*State_size, sizeof(double));
	if(S_middle==NULL || S_prime == NULL || S_dprime == NULL || S_e_prime == NULL || S_e_dprime == NULL || S_e_prime_s == NULL){
        #ifdef MATLAB_MEX_FILE
            mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
            fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }else{
        //provided there was memory for S_midle and S_prime, populate them
        for(i=0;i<S_length-1;i++){
            for(j=0;j<State_size;j++){
                S_middle[i*State_size+j]=(S[i*State_size+j]+S[(i+1)*State_size+j])/2;//S_middle = (S(:,2:end)+S(:,1:end-1))/2;
                S_prime[i*State_size+j] = (S[(i+1)*State_size+j]-S[i*State_size+j])/dtheta;//S_prime = (S(:,2:end)-S(:,1:end-1))/(dtheta);
				S_dprime[i*State_size + j] = 0;
            }
        }
		for (i = 0; i<S_length; i++){
			for (j = 0; j<State_size; j++){
				if (i < 2)
					S_e_prime[i*State_size + j] = (S[i*State_size + j] * (-11.0) / 6.0 + S[(i + 1)*State_size + j]*3.0 - S[(i + 2)*State_size + j] * 3.0 / 2.0 + S[(i + 3)*State_size + j] * 1.0 / 3.0) / dtheta;
				else if (i > (S_length - 3))
					S_e_prime[i*State_size + j] = (S[(i - 2)*State_size + j] * 1.0 / 2.0 - S[(i - 1)*State_size + j] * 2.0 + S[i*State_size + j] * 3.0 / 2.0) / dtheta;
				else
					S_e_prime[i*State_size + j] = (S[(i - 2)*State_size + j] * 1 / 12 - S[(i - 1)*State_size + j] * 2.0 / 3.0 + S[(i + 1)*State_size + j] * 2.0 / 3.0 - S[(i + 2)*State_size + j] * 1.0 / 12.0) / dtheta;
			}
		}
		for (i = 0; i<S_length; i++){
			for (j = 0; j<State_size; j++){
				if (i < 2)
					S_e_dprime[i*State_size + j] = (S[i*State_size + j] * (35.0) / 12.0 + S[(i + 1)*State_size + j] * -(26.0 / 3.0) + S[(i + 2)*State_size + j] * 19.0 / 2.0 - S[(i + 3)*State_size + j] * 14.0 / 3.0 + S[(i + 4)*State_size + j] * 11.0 / 12.0) / dtheta;
				else if (i >(S_length - 3))
					S_e_dprime[i*State_size + j] = (S[(i - 2)*State_size + j] * 1.0 - S[(i - 1)*State_size + j] * 2.0 + S[(i)*State_size + j] * 1.0) / dtheta;
				else
					S_e_dprime[i*State_size + j] = (S[(i - 2)*State_size + j] * (-1.0 / 12.0) + S[(i - 1)*State_size + j] * 4.0 / 3.0 - S[i*State_size + j] * 5.0 / 2.0 + S[(i + 1)*State_size + j] * 4.0 / 3.0 - S[(i + 2)*State_size + j] * 1.0 / 12.0) / dtheta;
			}
		}
		for (i = 0; i<S_length; i++){
			for (j = 0; j < State_size; j++){
				for (k = 0; k < State_size; k++){
					S_e_prime_s[i*State_size*State_size + j*State_size + k] = S_e_prime[i*State_size + j] * S_e_prime[i*State_size + k];
				}
			}
		}
    }
    
    //initialize b to speed_constant/norm(S_prime)
    if(*p_b == NULL)
        b = calloc(S_length,sizeof(double));
    else
        b  = realloc(*p_b, S_length*sizeof(double));
    a = calloc(S_length-1,sizeof(double));
    if(b==NULL || a == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }

    sum = 0;
    for(i=0;i<State_size;i++){
        sum += (S[State_size+i]-S[i])*(S[State_size+i]-S[i]);
    }
    if(p_params->initial_velocity < 0){
        #ifdef MATLAB_MEX_FILE
                mexPrintf("Negative initial velocity given, using the magnitude");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "Negative initial velocity given, using the magnitude\n");
        #endif
    }
    
    b_0 = p_params->initial_velocity*p_params->initial_velocity*dtheta*dtheta/sum;
    if(display){
        #ifdef MATLAB_MEX_FILE
                mexPrintf("b_0 is %f\n", b_0);
        #endif
        #ifndef MATLAB_MEX_FILE
                printf("b_0 is %f\n", b_0);
        #endif
    }
    b[0] = b_0;
    
    //initialize the rest of the b vector, provided it was allocated
    if(b != NULL){
        if(o_params == NULL || o_params->initial_b == NULL){
            for(i=1;i<S_length;i++){
                b[i] = b_0/2+dtheta*dtheta;
            }
        }else{
            memcpy(&b[1], &(o_params->initial_b[1]), (S_length-1)*sizeof(double));
        }
    }
    
    if(totime){
        diff = clock()-start;
        time = diff*1000.0/CLOCKS_PER_SEC;
        timers[0] = time;
        start = clock();
    }
    
    //set all of the pointers being passed to MakeA to NULL so that they will get allocated in MakeA
    csA21 = NULL;//a csparse storing (compressed) of the dynamics constraints matrix.
    indexA21_i = NULL;//the i indices of the sparse A21 array
    indexA21_j = NULL;//the j indices of the sparse A21 array
    indexA21_v = NULL;//the values of the sparse A21 array;
    b_Ax = NULL;
    
    G = NULL;//space for holding the gradient.
    Hc = NULL;//space for holding the Hessian.
    indexH_i = NULL;//space for storing the i indices of the Hessian
    indexH_j = NULL;//space for storing the j indices of the Hessian
    indexH_v = NULL;//space for storing teh values of the Hessian
    
    //set pointers to NULL so that you know what to free in the case of failure.
    Abig_c = NULL;
    Abig_i = NULL;
    Abig_j = NULL;
    dx = NULL;
    lu_S = NULL;
    
    if(status){
        A21index_length = so_MakeA(S, S_middle, S_prime, S_dprime, S_e_prime, S_e_dprime, S_e_prime_s, dtheta, U_size, S_length, State_size, &csA21, &indexA21_i, &indexA21_j, &indexA21_v, &b_Ax, b_0, variables, variables_length);
        //check that A21index_length executed properly
        if(A21index_length == -1)
            status = 0;
    }
        
    if(totime){//time required for A
       diff = clock()-start;
       time = diff*1000.0/CLOCKS_PER_SEC;
       timers[1] = time;
       start = clock();
    }
    
    //make x and u should occur after the make A21 so that I can make sure they are the correct size;
    x = calloc(csA21->n,sizeof(double));
    resid = calloc(csA21->m,sizeof(double));
    if(*p_u == NULL)
        u = calloc(csA21->n-2*(S_length-1),sizeof(double));
    else
        u  = realloc(*p_u, csA21->n-2*(S_length-1)*sizeof(double));
    if(*p_v == NULL)
        v = calloc(csA21->m,sizeof(double));
    else
        v = realloc(*p_v, csA21->m*sizeof(double));
    if(x == NULL || u == NULL || v == NULL || resid == NULL){
        #ifdef MATLAB_MEX_FILE
            mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
            fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }
    
    if(status){
        //check if an initial u has been provided, otherwise, initialize to zero
        if(o_params == NULL || o_params->initial_u == NULL){
            for(i=0;i<csA21->n-2*(S_length-1);i++){
                u[i] = 0;
            }
        }else{
            memcpy(u, o_params->initial_u, csA21->n-2*(S_length-1)*sizeof(double));
        }
        
        //initialize the v vector to zero
        for(i=0; i<csA21->m;i++)
            v[i] = 0;
        
        //initialize the b terms in x generated earlier.
        for(i=0;i<S_length-1;i++){
            x[i*(U_size+2)] = b[i+1];
            x[i*(U_size+2)+1] = (b[i+1]-b[i])/(2*dtheta);
            a[i] = x[i*(U_size+2)+1];
            for(j=0;j<U_size;j++)
                x[i*(U_size+2)+2+j] = u[i*U_size+j];//U initialization.
        }
        //for(i=0;i<csA21->n;i++){
        //mexPrintf("x[%d]:%f\n",i,x[i]);
        //}
        
        if(totime){//time required for more init
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[0] += time;
            start = clock();
        }
    }
    
    while(!solved && iterations < MAX_ITERATIONS && status){
        if(display){
            #ifdef MATLAB_MEX_FILE
                mexPrintf("on iteration number %d\n",iterations);
            #endif
            #ifndef MATLAB_MEX_FILE
                    printf("on iteration number %d\n",iterations);
            #endif
        }
        
        iterations += 1;
        
        //calculate the Hessian and Gradient
        Hindex_length = so_MakeHandG(b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, 0, &G, &Hc, &indexH_i, &indexH_j, &indexH_v,timers,totime,variables,variables_length);
        if(Hindex_length == -1)
            status = 0;
        else
            status = 1;
        /*for(i=1;i<Hindex_length;i++)
            mexPrintf("indexH_v[%d]=%e\n",i,indexH_v[i]);
        for(i=0;i<Hc->m;i++)
            mexPrintf("G[%d]=%e\n",i,G[i]);*/
        
        
        if(totime){//time for HandG
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[2] += time;
            start = clock();
        }
        
        //if the problem was infeasible last iteration, check if it feasible this iteration
        if(!feasible && status){
            //calculate A21*x-b_Ax;
            residual_new = 0;
            for(i=0;i<csA21->m;i++){
                resid[i] = -b_Ax[i];
                //mexPrintf("Resid[i]:%f\n",resid[i]);
            }
            cs_gaxpy(csA21, x, resid);
            for(i=0;i<csA21->m;i++){
                //mexPrintf("resid[i]:%f\n",resid[i]);
                residual_new += resid[i]*resid[i];
            }
            residual_new = sqrt(residual_new);
            
            if(display){
                #ifdef MATLAB_MEX_FILE
                        mexPrintf("Current Residual is %f\n",residual_new);
                #endif
                #ifndef MATLAB_MEX_FILE
                        printf("Current Residual is %f\n",residual_new);
                #endif
            }
            
            if(residual_new > residual){
                orderOkay = 0;
                if(display){
                    #ifdef MATLAB_MEX_FILE
                            mexPrintf("Sparsity structure changed");
                    #endif
                    #ifndef MATLAB_MEX_FILE
                            printf("Sparsity structure changed");
                    #endif
                }
            }
            
            residual = residual_new;
            
            if(residual < RESIDUAL_THRESHHOLD){
                feasible = 1;
                //note that this is the first time function evaluation is called, so this can be thought of like initializaiton.
                fx = so_function_evaluation(b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta,kappa, variables,variables_length);
                if(fx == -1){
                    status = 0;
                }
            }
            if(totime){//time required for feasibility check
                diff = clock()-start;
                time = diff*1000.0/CLOCKS_PER_SEC;
                timers[3] += time;
                start = clock();
            }
        }
        
        
        //Now generate the large A matrix for the algorithm step.
        //to do this, use the index arrays that are returned from make H and A12
        //the arrays are the same size on every iteration, so you only need to allocate them once;
        if(Abig_v == NULL && status){
            index_length = A21index_length*2+Hindex_length;
            Abig_i = calloc(index_length,sizeof(int));
            Abig_j = calloc(index_length,sizeof(int));
            Abig_v = calloc(index_length,sizeof(double));
            bbig = calloc(Hc->m+(S_length-1)*(State_size+1),sizeof(double));//bbig;
            dx = calloc((S_length-1)*(2+U_size),sizeof(double));//dx from the solution;
            if(Abig_i == NULL || Abig_j == NULL || Abig_v == NULL || bbig == NULL || dx == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("insufficient memory");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "insufficient memory\n");
                #endif
                status = 0;
            }
        }
        
        //stuff the matrix Abig.
        if(status){
            so_stuffing(Abig_v, Abig_i, Abig_j, indexH_v, indexH_i, indexH_j, Hindex_length, indexA21_v, indexA21_i, indexA21_j, A21index_length, Hc->m, S_length, State_size, U_size);
            
            
            //construct the large sparse array
            Abig.nzmax = index_length;
            Abig.i = Abig_i;
            Abig.p = Abig_j;
            Abig.x = Abig_v;
            Abig.nz = Abig.nzmax;
            Abig.m = Hc->m+csA21->m;
            Abig.n = Hc->n+csA21->m;
            if(totime){
                start2 = clock();
            }
            Abig_c = cs_triplet(&Abig);
            if(totime){
                diff = clock()-start2;
                time = diff*1000.0/CLOCKS_PER_SEC;
                timers[18] += time;
            }
            if(Abig_c == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("matrix compression failed, out of memory");
                #endif
                        #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "matrix compression failed, out of memory\n");
                #endif
                status = 0;
            }
        }

        if(status){
            cs_dupl(Abig_c);
            cs_dropzeros(Abig_c);
            
            //finish creating b;
            for(i=0;i<(2+U_size)*(S_length-1);i++){
                bbig[i] = -G[i];//-G potion
            }
            
            if(!feasible){//if not feasible the rest of this array should be zeros and will be set earlier.
                for(i=0;i<csA21->m;i++)
                    bbig[i+(2+U_size)*(S_length-1)]= -resid[i];
            }else{
                //since the solve overwrites bbig, we now need to set the rest of the vector to zero each time if it is zero;
                for(i=0;i<csA21->m;i++)
                    bbig[i+(2+U_size)*(S_length-1)]= 0;
            }
            
            if(totime){//time required to create the large matrix
                diff = clock()-start;
                time = diff*1000.0/CLOCKS_PER_SEC;
                timers[4] += time;
                start = clock();
            }           			

            orderOkay = so_lu_solve(Abig_c, bbig, &lu_S, orderOkay);
            if(orderOkay==0){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("solve step failed");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "solve step failed\n");
                #endif
                status = 0;
            }
            
            if(totime){//time for solve
                diff = clock()-start;
                time = diff*1000.0/CLOCKS_PER_SEC;
                timers[5] += time;
                start = clock();
            }
        }
        
        if(status){
            //copy data into the dx array;
            memcpy(dx, bbig, (S_length-1)*(2+U_size)*sizeof(double));
            if(!feasible){
                status = so_linesearch_infeasible(&bbig[(S_length-1)*(2+U_size)], v, x, dx, b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, beta, alpha, G, csA21, b_Ax, variables, variables_length);
                if(totime){//time for infeasible line search
                    diff = clock()-start;
                    time = diff*1000.0/CLOCKS_PER_SEC;
                    timers[6] += time;
                    start = clock();
                }
            }else{
                
                //calculate Lambda = dx'*H*dx;
                data = calloc((S_length-1)*(2+U_size), sizeof(double));
                if(data == NULL){
                    #ifdef MATLAB_MEX_FILE
                            mexErrMsgTxt("insufficient memory");
                    #endif
                            #ifndef MATLAB_MEX_FILE
                            fprintf(stderr, "insufficient memory\n");
                    #endif
                    status = 0;
                }
                if(status){
                    //initialize data to be zero;
                    for(i=0;i<(S_length-1)*(2+U_size);i++)
                        data[i]=0;
                    cs_gaxpy(Hc, dx, data);
                    Lambda = 0;
                    for(i=0;i<(S_length-1)*(2+U_size);i++){
                        Lambda += data[i]*dx[i];
                        //mexPrintf("dx[%d]=%e\n",i,dx[i]);
                    }
                    free(data);
                    if(display){
                        #ifdef MATLAB_MEX_FILE
                                mexPrintf("Lambda is %f\n", Lambda);
                        #endif
                                #ifndef MATLAB_MEX_FILE
                                printf("Lambda is %f\n", Lambda);
                        #endif
                    }
                }
                
                if(status){
                    //check Lambda to see if you're done and if kappa needs to be reduced.
                    if(Lambda/2 < epsilon){
                        if(kappa*S_length > epsilon && (flags_kappa == 1)){
                            kappa = kappa/mu;
                            fx = so_function_evaluation(b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, variables, variables_length);
                            if(fx == -1)
                                status = 0;
                            if(display){
                                #ifdef MATLAB_MEX_FILE
                                        mexPrintf("reducing kappa\n");
                                #endif
                                        #ifndef MATLAB_MEX_FILE
                                        printf("reducing kappa\n");
                                #endif
                            }
                        }else{
                            solved = 1;
                            //mexPrintf("The problem is solved!\n");
                        }
                    }else{
                        //perform the feasible line search
                        fx = so_linesearch_feasible(&bbig[(S_length-1)*(2+U_size)], v, x, dx, b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, G, alpha, beta, fx, variables, variables_length);
                        if(fx == -1)
                            status = 0;
                    }
                }
                if(totime){//time for feasible line search
                    diff = clock()-start;
                    time = diff*1000.0/CLOCKS_PER_SEC;
                    timers[7] += time;
                    start = clock();
                }
            }
        }
        
        //free the large A matrix you create
        if(Abig_c != NULL){
            cs_spfree(Abig_c);
            Abig_c = NULL;
        }
        
        if(totime){//time counted as part of making A
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[1] += time;
            start = clock();
        }
    }

    //Set all of the values to be returned
    *p_b = b;
    *p_u = u;
    *p_v = v;
    *p_fx = fx;
    *p_iterations = iterations;
    
    if(S_middle != NULL)
        free(S_middle);
    if(S_prime!= NULL)
        free(S_prime);
    if(S_dprime != NULL)
        free(S_dprime);
	if (S_e_dprime != NULL)
		free(S_e_dprime);
	if (S_e_prime != NULL)
		free(S_e_prime);
	if (S_e_prime_s != NULL)
		free(S_e_prime_s);


    //free all of the space allocated for csA21;
    if(indexA21_i != NULL)
        free(indexA21_i);
    if(indexA21_j != NULL)
        free(indexA21_j);
    if(indexA21_v != NULL)
        free(indexA21_v);
    if(b_Ax != NULL)
        free(b_Ax);
    if(csA21 != NULL)
        cs_spfree(csA21);
    
    if(a != NULL)
        free(a);
    if(x !=  NULL)
        free(x);
    //destory arrays that are constant for HG
    if(G != NULL)
        free(G);
    
    if(Hc != NULL)
        cs_spfree(Hc);
    if(indexH_i != NULL)
        free(indexH_i);
    if(indexH_j != NULL)
        free(indexH_j);
    if(indexH_v != NULL)
        free(indexH_v);
    
    //clear the elements created for the use in calculating the residual;
    if(resid != NULL)
        free(resid);
    
    //free space created to make Abig
    if(Abig_i !=NULL)
        free(Abig_i);
    if(Abig_j != NULL)
        free(Abig_j);
    if(Abig_v != NULL)
        free(Abig_v);
    if(bbig != NULL)
        free(bbig);
    if(dx != NULL)
        free(dx);
    if(lu_S != NULL)
        cs_sfree(lu_S);
    
    if(totime){//time for solve
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[8] = time;
    }
    
   return status; 
}

int so_calcHbarrier(double* H_barrier, int S_length, int U_size, int** pindex_i, int** pindex_j, double** pindex_v){
    /*so_calcHBarrier, calculates the portion of the Hessian produced by the user defined barrier function.
 *BarrierArray: the three dimensional array of inputs produced by the barrier function
 *index_i: an array of integers to store the i indices
 *index_j: an array of integers to store the j indices
 *index_v: an array of doubles to store the values corresponding to those indices;
 *return: the length of the index arrays, or -1 if failure;
*/
    int i, n, m, block_length, index_length;
    int status;
    double *index_v;
    int *index_i, *index_j;
    
    status = 1;
    
    block_length = (2+U_size)*(2+U_size);//the size of each H block.
    index_length = block_length*(S_length-1);//the length of H_Barrier.
    
    //note that to speed things up further, the i and j calculations need only be performed once, and could be stored in a different implementation.
    if(*pindex_i == NULL)
        index_i = calloc(index_length,sizeof(int));
    else
        index_i = realloc(*pindex_i,index_length*sizeof(int));
    
    if(*pindex_j == NULL)
        index_j = calloc(index_length,sizeof(int));
    else
        index_j = realloc(*pindex_j,index_length*sizeof(int));
    
    if(*pindex_v == NULL)
        index_v = calloc(index_length,sizeof(double));
    else
        index_v = realloc(*pindex_v,index_length*sizeof(double));
    
    if(index_i == NULL || index_j == NULL || index_v == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
                status = 0;
    }
    *pindex_i = index_i;
    *pindex_j = index_j;
    *pindex_v = index_v;
    
    for(i=0;i<S_length-1;i++){
        //copy the H from one time step into TempH so that it can be manipulated to account for the averaging of b across timesteps.
        for(n=0;n<(U_size+2);n++){
            //populate the indices of the block matrix.
            for(m=0;m<(U_size+2);m++){
                index_i[i*block_length+n*(U_size+2)+m]=m+i*(U_size+2);
                index_j[i*block_length+n*(U_size+2)+m]=n+i*(U_size+2);
            }
        }
    }
    //copy the derivative data into index_v;
    memcpy(index_v,H_barrier,index_length*sizeof(double));
    if(status)
        return(index_length);
    else
        return -1;
}

int so_MakeHandG(double* b, double* S_middle, double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, int flag_calculate, double** G, cs** pH, int** pindexH_i, int** pindexH_j, double** pindexH_v, double* timers, int istimed, double* variables, int variables_length){
    /*b: an array of the b values
     *S_middle: the S values evaluated at their midpoints
     *u: an array of the u values as determined at the midpoints of the S values
     *dtheta: the distance along the path between steps
     *kappa: the weighting term for the barrier
     *barrier: the function handle used to call the barrier function
     *flag_calculate: a flag showing which terms (H and G) need to be calculated: 0 for both 1 for H, 2 for G
     *G: a pointer to an array for storing the gradient
     *H: a pointer to a csparse structure in which to store the hessian
     *indexH_i: a pointer to an array of ints to store the i index to make the KKT system with
     *indexH_j: a pointer to an array of ints store the j index to make the KKT system with
     *indexH_v: a pointer to an array of doubles to store the values corresponding with the indices in indexH_i and indexH_j for the creation of the KKT system
     *return:  the length of the indexH arrays.*/
    double time;
    int doG, doH, status, statusB, statusF, length_Hb;
    int i;//index for loops.
    int *indexH_i, *indexH_j, index_length, *indexHb_i, *indexHb_j;
    double *data, *bsqrt, *bsqrt_mid, *bsqrt_mid2, *bsqrt3, *bsqrt_mid3, *indexH_v, *indexHb_v, *H_barrier, *G_barrier, *F_barrier;
    cs H, *Hc;//Hc is compressed form H.
    clock_t start, diff;
    if(istimed){
        start = clock();
        if(timers == NULL)
            istimed = 0;
    }
    
    status = 1;//checks memory status
    statusB = 1;
    statusF = 0;
    index_length = 3*S_length-5+(U_size+2)*(U_size+2)*(S_length-1);
    H_barrier = NULL;
    G_barrier = NULL;
    F_barrier = NULL;
    indexHb_i = NULL;
    indexHb_j = NULL;
    indexHb_v = NULL;
    //set to NULL so we know if they need to be freed in the case of memory failure
    bsqrt = NULL;
    bsqrt_mid = NULL;
    bsqrt_mid2 = NULL;
    
    //check to make sure that b>0 if not, then terminate
    for(i=1;statusB && i < S_length;i++){
        if(b[i]<0)
            statusB = 0;
    }
    
    if(statusB){
        doG = (flag_calculate == 0) || (flag_calculate == 2);
        doH = (flag_calculate == 0) || (flag_calculate == 1);
    }else{
        doG = 0;
        doH = 0;
    }

    if(istimed){
        diff = clock()-start;
        time = diff*1000.0/CLOCKS_PER_SEC;
        timers[9] += time;//HG initialization;
        start = clock();
    }
    
    if(!doG){
        if(*G == NULL){
            *G = calloc(1, sizeof(double));
            if(*G == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("insufficient memory");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "insufficient memory\n");
                #endif
                status = 0;
            }
        }
        if(status)
            *G[0] = INFINITY;
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[10] += time;//G infeasible;
            start = clock();
        }
    }
    
    if(statusB){
        bsqrt = malloc(S_length*sizeof(double));
        bsqrt_mid = malloc((S_length-1)*sizeof(double));
        bsqrt_mid2 = malloc((S_length-1)*sizeof(double));
        if(bsqrt == NULL || bsqrt_mid == NULL || bsqrt_mid2 == NULL){
            #ifdef MATLAB_MEX_FILE
                    mexErrMsgTxt("insufficient memory");
            #endif
                    #ifndef MATLAB_MEX_FILE
                    fprintf(stderr, "insufficient memory\n");
            #endif
                    status = 0;
        }else{
            //calculates some values used by everyone.
            for(i=0;i<S_length;i++){
                bsqrt[i] = sqrt(b[i]);
                if(i!=0){
                    bsqrt_mid[i-1] = bsqrt[i-1]+bsqrt[i];
                    bsqrt_mid2[i-1] = 1/(bsqrt_mid[i-1]*bsqrt_mid[i-1]);
                }
            }
        }
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[11] += time;//bqrt generation;
            start = clock();
        }
    }
    
    //calculates the terms in the Hessian and Gradient dependent on problem formulation, not on the barrier.  This may need to be decremented by 1.
    if(doH && status){
        bsqrt3 = malloc((S_length-1)*sizeof(double));
        bsqrt_mid3 = malloc((S_length-1)*sizeof(double));
        //3*S_length-5 is for this term, the remainder of the space allocated is for the ouput from CalcH.
        if(*pindexH_i==NULL)
            indexH_i = malloc(index_length*sizeof(int));
        else
            indexH_i = realloc(*pindexH_i, index_length*sizeof(int));//these reallocs may not be necessary, since they should be the same size.
        
        if(*pindexH_j==NULL)
            indexH_j = malloc(index_length*sizeof(int));
        else
            indexH_j = realloc(*pindexH_j, index_length*sizeof(int));
        
        if(*pindexH_v==NULL)
            indexH_v = malloc(index_length*sizeof(double));
        else
            indexH_v = realloc(*pindexH_v, index_length*sizeof(double));
        
        if(bsqrt3 == NULL || bsqrt_mid3 == NULL || indexH_i == NULL || indexH_j == NULL || indexH_v == NULL){
            #ifdef MATLAB_MEX_FILE
                    mexErrMsgTxt("insufficient memory");
            #endif
            #ifndef MATLAB_MEX_FILE
                    fprintf(stderr, "insufficient memory\n");
            #endif
            status = 0;
        }else{
            //generate the Hessian necessary for the minimum time portion of the problem.
            
            
            for(i=0;i<S_length-1;i++){
                //generate some terms necessary for the calculation of the Hessian for the minimum time component.
                bsqrt3[i] = 1/(bsqrt[i+1]*bsqrt[i+1]*bsqrt[i+1]);
                bsqrt_mid3[i] = 1/(bsqrt_mid[i]*bsqrt_mid[i]*bsqrt_mid[i]);
                //generate the diagonal elements of the Hessian.
                indexH_v[i] = dtheta*(bsqrt3[i]*bsqrt_mid2[i]/2+bsqrt_mid3[i]/b[i+1])+kappa/(b[i+1]*b[i+1]);
                indexH_i[i] = (2+U_size)*i;
                indexH_j[i] = (2+U_size)*i;
                //generate the offdiagonal elements of the Hessian, that arise because we are looking at the midpoint between bs.
                if(i>0){//for all except the end term, so all elements are i-1, so that the b terms can be calculated in the same loop.
                    indexH_v[i-1]+= dtheta*(bsqrt3[i-1]*bsqrt_mid2[i-1]/2+1/b[i]*bsqrt_mid3[i-1]);
                    indexH_v[i+S_length-2] = dtheta*(bsqrt_mid3[i]/(bsqrt[i]*bsqrt[i+1]));
                    indexH_i[i+S_length-2] = (2+U_size)*(i-1);
                    indexH_j[i+S_length-2] = (2+U_size)*(i);
                    indexH_v[i+2*S_length-4] = dtheta*(bsqrt_mid3[i]/(bsqrt[i]*bsqrt[i+1]));
                    indexH_i[i+2*S_length-4] = (2+U_size)*(i);
                    indexH_j[i+2*S_length-4] = (2+U_size)*(i-1);
                }
            }
        }
        if(bsqrt3 != NULL)
            free(bsqrt3);
        if(bsqrt_mid3 != NULL)
            free(bsqrt_mid3);
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[12] += time;//H problem formulation;
            start = clock();
        }
    }
    
    
    //perform setup necessary to call the barrier function
    if(statusB && status){//if the B status already failed then you aren't going to be doing anything.
        if(flag_calculate<2){
            H_barrier = malloc((2+U_size)*(2+U_size)*(S_length-1)*sizeof(double));
            if(H_barrier == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("Insufficient memory to allocate the hessian for the barrier");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "Insufficient memory to allocate the hessian for the barrier\n");
                #endif
                status = 0;
            }
        }
        if(flag_calculate == 0 || flag_calculate == 2){
            G_barrier = malloc((2+U_size)*(S_length-1)*sizeof(double));
            if(G_barrier == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("Insufficient memory to allocate the gradient for the barrier");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "Insufficient memory to allocate the gradient for the barrier\n");
                #endif
                status = 0;
            }
        }
        if(status){
            statusF = so_barrier(S_middle,S_prime, S_dprime, &b[1],a,u, S_length,U_size,State_size, flag_calculate,kappa,variables,variables_length,H_barrier,G_barrier,F_barrier);
            if(statusF == -1){
                status = 0;
                statusF = 0;
            }
        }
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[13] += time;//barrier function;
            start = clock();
        }
    }
    
    
    if(doG && status){
        if(statusF == 0){
            if(*G == NULL){
                *G = calloc(1, sizeof(double));
                if(*G == NULL){
                    #ifdef MATLAB_MEX_FILE
                            mexErrMsgTxt("Insufficient memory");
                    #endif
                    #ifndef MATLAB_MEX_FILE
                            fprintf(stderr, "Insufficient memory\n");
                    #endif
                    status = 0;
                }
            }
            *G[0] = INFINITY;
        }else{
            //allocate space for G if it has not yet been allocated.
            if(*G == NULL){
                *G = calloc((2+U_size)*(S_length-1), sizeof(double));
                if(*G == NULL){
                    #ifdef MATLAB_MEX_FILE
                            mexErrMsgTxt("Insufficient memory");
                    #endif
                    #ifndef MATLAB_MEX_FILE
                            fprintf(stderr, "Insufficient memory\n");
                    #endif
                    status = 0;
                }
                data = *G;
            }else{
                *G = realloc(*G, (2+U_size)*(S_length-1)*sizeof(double));
                if(*G == NULL){
                    #ifdef MATLAB_MEX_FILE
                            mexErrMsgTxt("Insufficient memory");
                    #endif
                    #ifndef MATLAB_MEX_FILE
                            fprintf(stderr, "Insufficient memory\n");
                    #endif
                    status = 0;
                }
                data = *G;
            }
            if(status){
                memcpy(data, G_barrier, (2+U_size)*(S_length-1)*sizeof(double));
                for(i=0;i<S_length-1;i++){
                    //populate the G1 terms which correspond to the minimum time calculation.
                    data[i*(2+U_size)] += dtheta*(-bsqrt_mid2[i]/bsqrt[i+1])-kappa/b[i+1];
                    if(i!=S_length-2){
                        data[i*(2+U_size)] += -dtheta*bsqrt_mid2[i+1]/bsqrt[i+1];
                    }
                }
            }
        }
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[14] += time;//G formulation;
            start = clock();
        }
    }
    
    
    if(doH && status){
        length_Hb = so_calcHbarrier(H_barrier, S_length, U_size, &indexHb_i, &indexHb_j, &indexHb_v);
        if(length_Hb == -1)
            status = 0;
        
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[15] += time;//calcH;
            start = clock();
        }
       
        if(status){
            memcpy(&indexH_v[3*S_length-5], indexHb_v, length_Hb*sizeof(double));
            memcpy(&indexH_i[3*S_length-5], indexHb_i, length_Hb*sizeof(int));
            memcpy(&indexH_j[3*S_length-5], indexHb_j, length_Hb*sizeof(int));
            
            H.nzmax = 3*S_length-5+length_Hb;
            H.i = indexH_i;
            H.p = indexH_j;
            H.x = indexH_v;
            H.nz = H.nzmax;
            H.m = (S_length-1)*(2+U_size);
            H.n = (S_length-1)*(2+U_size);
            Hc = cs_triplet(&H);
            if(Hc == NULL){
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("matrix compression failed, out of memory");
                #endif
                        #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "matrix compression failed, out of memory\n");
                #endif
                        status = 0;
            }
            cs_dupl(Hc);
            cs_dropzeros(Hc);
            //Now I just need away to return this matrix.
            if(*pH != NULL)
                cs_spfree(*pH);
            
            *pH = Hc;
            
            //make sure that all of the index arrays are returned
            *pindexH_v = indexH_v;
            *pindexH_i = indexH_i;
            *pindexH_j = indexH_j;
        }
        if(indexHb_i != NULL)
            free(indexHb_i);
        if(indexHb_j != NULL)
            free(indexHb_j);
        if(indexHb_v != NULL)
            free(indexHb_v);
        if(istimed){
            diff = clock()-start;
            time = diff*1000.0/CLOCKS_PER_SEC;
            timers[16] += time;//big H generation;
            start = clock();
        }
    }
    
    if(statusB){
        if(H_barrier!=NULL)
            free(H_barrier);
        if(G_barrier!=NULL)
            free(G_barrier);
        if(F_barrier!=NULL)
            free(F_barrier);
        if(bsqrt != NULL)
            free(bsqrt);
        if(bsqrt_mid != NULL)
            free(bsqrt_mid);
        if(bsqrt_mid2 != NULL)
            free(bsqrt_mid2);
    }
    
    
    
    if(istimed){
        diff = clock()-start;
        time = diff*1000.0/CLOCKS_PER_SEC;
        timers[17] += time;//freeing;
    }
    if(status)
        return index_length;
    else
        return -1;
}

int so_linesearch_infeasible(double* w, double* v, double* x, double* dx, double* b, double* S_middle, double* S_prime, double* S_dprime, double *u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double beta, double alpha, double* G, cs* A21, double* b_Ax, double* variables, int variables_length) {
    //inputs are w [0],v [1],x [2],dx [3], S_middle [4], dtheta [5], kappa [6], beta [7], G [8], A_21 [9], b_Ax[10], barrier [11], alpha[12];
    //plhs[0] = X, plhs[1] = V, plhs[2] = b, because it is needed later;
    double *X, *V, *term1, *term2, *Gnew;
    double t, LHS, RHS;
    int MAX_ITERATIONS, iterations;
    int i, j, status, doloop, calcRHS;
    int* indexH_i, *indexH_j;
    double *indexH_v;
    cs *A21T_c, *Htrash;
    
    t = 1;
    doloop = 1;
    calcRHS = 1;
    Gnew = NULL;
    Htrash = NULL;
    indexH_i = NULL;
    indexH_j = NULL;
    indexH_v = NULL;
    A21T_c = NULL;
    term1 = NULL;
    term2 = NULL;
    MAX_ITERATIONS = 100;
    iterations = 0;
    status = 1;
    
    X = calloc(A21->n,sizeof(double));
    V = calloc(A21->m,sizeof(double));
    if(X == NULL || V == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }
    
    if(status){
        do{
            //compute X = x+t*dx;
            for(i=0;i<A21->n;i++){
                X[i] = x[i]+t*dx[i];
            }
            
            //compute b and u (or just copy from the computation);
            for(i=0;i<S_length-1;i++){
                b[i+1] = X[i*(U_size+2)];
                a[i] = X[i*(U_size+2)+1];
                for(j=0;j<U_size;j++){
                    u[i*U_size+j]=X[i*(U_size+2)+2+j];
                }
            }
            
            //calculate V
            for(i=0;i<A21->m;i++){
                V[i] = v[i]*(1-t)+t*w[i];
            }
            
            status = so_MakeHandG(b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, 2, &Gnew, &Htrash, &indexH_i, &indexH_j, &indexH_v, NULL, 0, variables, variables_length);//the last three are H_i, H_j, H_v
            if(status == -1)
                status = 0;
            else
                status = 1;
            
            if(status && !(G[0]==INFINITY)){
                //perform right hand calculation, but only on the first time through;
                if(calcRHS){
                    RHS = 0;
                    //calculate A21*x-b_Ax;
                    term1 = calloc(A21->m, sizeof(double));
                    if(term1 == NULL){
                        #ifdef MATLAB_MEX_FILE
                                mexErrMsgTxt("insufficient memory");
                        #endif
                        #ifndef MATLAB_MEX_FILE
                                fprintf(stderr, "insufficient memory\n");
                        #endif
                        status = 0;
                    }
                    if(status){
                        for(i=0;i<A21->m;i++)
                            term1[i] = -b_Ax[i];
                        cs_gaxpy(A21, x, term1);
                        for(i=0;i<A21->m;i++)
                            RHS += term1[i]*term1[i];
                        
                        //calculate G+A21'*v
                        A21T_c = cs_transpose(A21, 1);
                        if(A21T_c == NULL){
                            #ifdef MATLAB_MEX_FILE
                                    mexErrMsgTxt("Matrix transpose failed, out of memory");
                            #endif
                            #ifndef MATLAB_MEX_FILE
                                    fprintf(stderr, "Matrix transpose failed, out of memory\n");
                            #endif
                            status = 0;
                        }
                        term2 = calloc(A21->n, sizeof(double));
                        if(term2 == NULL){
                            #ifdef MATLAB_MEX_FILE
                                    mexErrMsgTxt("insufficient memory");
                            #endif
                            #ifndef MATLAB_MEX_FILE
                                    fprintf(stderr, "insufficient memory\n");
                            #endif
                            status = 0;
                        }
                    }
                    if(status){
                        memcpy(term2, G, A21->n*sizeof(double));
                        cs_gaxpy(A21T_c, v, term2);
                        for(i=0;i<A21->n;i++)
                            RHS += term2[i]*term2[i];
                        calcRHS = 0;
                    }
                }
                LHS = 0;
                
                if(status && !(Gnew[0]==INFINITY)){
                    //Calculate A21*X-b_Ax;
                    for(i=0;i<A21->m;i++)
                        term1[i] = -b_Ax[i];
                    cs_gaxpy(A21, X, term1);
                    for(i=0;i<A21->m;i++)
                        LHS += term1[i]*term1[i];
                    
                    //Gnew+A21'*V
                    memcpy(term2, Gnew, A21->n*sizeof(double));
                    cs_gaxpy(A21T_c, V, term2);
                    for(i=0;i<A21->n;i++)
                        LHS+= term2[i]*term2[i];
                    
                    //perform left hand calculation.
                    //mexPrintf("%e<%e\n",LHS,(1-alpha*t)*(1-alpha*t)*RHS);
                    if(LHS<(1-alpha*t)*(1-alpha*t)*RHS)
                        doloop = 0;
                }
            }
            
            if(doloop && iterations > MAX_ITERATIONS){
                doloop = 0;
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("The backward linesearch loop did not terminate correctly");
                #endif
                #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "The backward linesearch loop did not terminate correctly\n");
                #endif
                status = 0;
            }
            iterations += 1;
            t *= beta;
        }while(doloop && status);
    }

    //mexPrintf("Gnew[0] is %f\n",Gnew[0]);
    //copy X and V into x and v before freeing;
    if(status){
        memcpy(x, X, A21->n*sizeof(double));
        memcpy(v, V, A21->m*sizeof(double));
    }
    if(X != NULL)
        free(X);
    if(V != NULL)
        free(V);
    if(Gnew != NULL)
        free(Gnew);
    if(A21T_c != NULL)
        cs_spfree(A21T_c);
    if(term1 != NULL)
        free(term1);
    if(term2 != NULL)
        free(term2);
    if(indexH_i != NULL)
        free(indexH_i);
    if(indexH_j != NULL)
        free(indexH_j);
    if(indexH_v != NULL)
        free(indexH_v);
    if(Htrash != NULL)
        cs_spfree(Htrash);
   
    return status;
}

double so_function_evaluation(double* b, double* S_middle, double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double* variables, int variables_length){
    int i, b_okay, status, valid;
    double *bsqrt, *bsqrt_mid, *F;
    double fx;
    
    fx = 0;
    b_okay = 1;
    bsqrt = calloc(S_length,sizeof(double));
    bsqrt_mid = calloc(S_length-1,sizeof(double));
    status = 1;
    if(bsqrt == NULL || bsqrt_mid == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }
    
    for(i=0;i<S_length && b_okay;i++){
        if(b[i]<0)//if b[i] < zero trajectory is infeasible and value should be inf.
            b_okay = 0;
        else{//calculate bsqrt and bsqrtmid for use in the funciton evaluation.  (in particular the b>0 barrier)
            bsqrt[i] = sqrt(b[i]);
            if(i!=0)
                bsqrt_mid[i-1] = bsqrt[i-1]+bsqrt[i];
        }
    }
    
    if(b_okay && status){
        //call the barrier function
        //for user easer, the memory is allocated here so the user doesn't need to worry about memory allocation,
        //However, this memory really only needs to be allocated once, so perhaps move up to an even higher level.
        F = malloc((S_length-1)*sizeof(double));
        valid = so_barrier(S_middle,S_prime, S_dprime, &b[1],a,u, S_length,U_size,State_size, 3,kappa,variables,variables_length,NULL,NULL,F);
        if(valid == -1)
            status = 0;
        else
            status = 1;

        if(valid && status){
            for(i=0;i<S_length-1;i++){
                fx += 2.0*dtheta/bsqrt_mid[i]+F[i]-kappa*(log(b[i+1]));
            }
        }else
            b_okay = 0;
        //free the results of the barrier function
        if(F!=NULL)
            free(F);
    }
    
    //free the sqrt arrays that were just created.
    if(bsqrt != NULL)
        free(bsqrt);
    if(bsqrt_mid != NULL)
        free(bsqrt_mid);
    
    if(status){
        if(b_okay)
            return(fx);
        else
            return(INFINITY);
    }else
        return -1;
}

double so_linesearch_feasible(double* w, double* v, double* x, double* dx, double* b, double* S_middle, double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double* G, double alpha, double beta, double fx, double* variables, int variables_length){
    double t, RHS, fX;
    double* X;
    int status, i, j, iterations, ITERATIONS_MAX, doloop, calcRHS;
    
    status = 1;
    t = 1;
    iterations = 0;
    ITERATIONS_MAX = 100;
    doloop = 1;
    calcRHS = 1;
    
    X = calloc((S_length-1)*(2+U_size),sizeof(double));
    if(X == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
        status = 0;
    }
    
    if(status){
        do{
            //compute X = x+t*dx;
            for(i=0;i<(S_length-1)*(2+U_size);i++){
                X[i] = x[i]+t*dx[i];
            }
            //compute b, u (or just copy from the computation);
            for(i=0;i<S_length-1;i++){
                b[i+1] = X[i*(U_size+2)];
                a[i] = X[i*(U_size+2)+1];
                for(j=0;j<U_size;j++){
                    u[i*U_size+j]=X[i*(U_size+2)+2+j];
                }
            }
            
            fX = so_function_evaluation(b, S_middle, S_prime, S_dprime, u, a, S_length, U_size, State_size, dtheta, kappa, variables, variables_length);
            if(fX == -1){//error occurred
                status = 0;
            }
            
            if(calcRHS){
                RHS = 0;
                //mexPrintf("G[0]=%f,dx[0]=%f\n",data[0],data2[0]);
                for(i=0;i<(S_length-1)*(U_size+2);i++){
                    RHS += G[i]*dx[i];
                }
                //mexPrintf("\n");
                RHS *= alpha;
                calcRHS = 0;
            }
            
            //mexPrintf("%e<%e\n",fX,fx+t*RHS);
            if(fX != INFINITY && fX < fx+t*RHS)
                doloop = 0;
            else{
                iterations += 1;
                t *= beta;
            }
            
            if(iterations > ITERATIONS_MAX){
                doloop = 0;
                #ifdef MATLAB_MEX_FILE
                        mexErrMsgTxt("Feasible Linesearch failed to exit loop");
                #endif
                        #ifndef MATLAB_MEX_FILE
                        fprintf(stderr, "Feasible Linesearch failed to exit loop\n");
                #endif
                status = 0;
            }
            
        }while(doloop && status);
    }
    
    //update X
    memcpy(x,X,(S_length-1)*(U_size+2)*sizeof(double));
    if(X != NULL)
        free(X);
    
    if(status){
        //update V
        for(i=0;i<(S_length-1)*(State_size+1);i++){
            v[i] = v[i]*(1-t)+t*w[i];
        }
        
        
        return(fX);
    }else
        return -1;
}

int so_MakeA(double* S, double* S_middle, double* S_prime, double* S_dprime, double* S_e_prime, double* S_e_dprime, double* S_e_prime_s, double dtheta, int U_size, int S_length, int State_size, cs** pAc, int** pA_i, int** pA_j, double** pA_v, double** pb_Ax, double b_0, double* variables, int variables_length){
    //generates the matrices necessary to represent the dynamics
    //if things fail, return -;
    double *A_v, *b_Ax;
    double *datac1, *datac2, *datam;
    double *M_dynamics, *R_dynamics, *C_dynamics, *d_dynamics;
    volatile int *A_i, *A_j, *indexR_i, *indexR_j;
    int block_length, index_length;//size parameters
    int i, j, status;//indices for loops
    cs A;
	ptrdiff_t m_m, n_m, p_m, m_m2, n_m2;//variables for the use of mtimes;
    char *chn = "N";//our matrices for our matrix multiplies are in column order, so no transpose necessary.
    double one = 1.0, zero = 0.0;//values for use in the matrix muliply one*A*B+zero*C
	FILE *fd;

    //a flag to make sure everything has executed properly
    status = 1;
    
    m_m = State_size;
    n_m = State_size;
    p_m = 1;
	m_m2 = State_size;
	n_m2 = State_size*State_size;

    //generate the first and second derivative of S with respect to theta necessary for the algorithm.
    //S_prime = calloc((S_length-1)*State_size, sizeof(double));
    R_dynamics = malloc(S_length*U_size*State_size*sizeof(double));
    M_dynamics = malloc(S_length*State_size*State_size*sizeof(double));
	C_dynamics = malloc(S_length*State_size*State_size*State_size*sizeof(double));
    d_dynamics = malloc(S_length*State_size*sizeof(double));
    indexR_i = calloc(U_size*State_size,sizeof(int));
    indexR_j = calloc(U_size*State_size,sizeof(int));
    datam  = calloc(m_m*n_m,sizeof(double));
    datac1 = calloc(m_m*n_m,sizeof(double));
    datac2 = calloc(m_m*n_m,sizeof(double));
    if(R_dynamics == NULL || M_dynamics == NULL || C_dynamics == NULL || d_dynamics == NULL || indexR_i == NULL || indexR_j == NULL ||datam == NULL || datac1 == NULL || datac2 == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("Insufficient Memory for allocation of dynamics data");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "Insufficient Memory for allocation of dynamics data\n");
        #endif
        status = 0;
    }
    
    //allocate space for the A matrix;
    block_length = (State_size*3+State_size*U_size+2);//4 for a,c*2;R;b'=a;
    index_length = block_length*(S_length-1)-State_size;
    if(*pA_i == NULL)
        A_i = calloc(index_length, sizeof(int));//this is written starting from rest...and thus b(0) is 1.
    else
        A_i = realloc(*pA_i, index_length*sizeof(int));
    
    if(*pA_j == NULL)
        A_j = calloc(index_length, sizeof(int));
    else
        A_j = realloc(*pA_j,index_length*sizeof(int));
    
    if(*pA_v == NULL)
        A_v = calloc(index_length, sizeof(double));
    else
        A_v = realloc(*pA_v,index_length*sizeof(double));
    
    if(*pb_Ax == NULL)
        b_Ax = calloc((S_length-1)*(State_size+1),sizeof(double));
    else
        b_Ax = realloc(*pb_Ax,(S_length-1)*(State_size+1)*sizeof(double));
    
    if(A_i == NULL || A_j == NULL || A_v == NULL || b_Ax == NULL){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("insufficient memory");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "insufficient memory\n");
        #endif
                status = 0;
    }
    
    
    
    status = so_dynamics(S, S_prime, S_length, State_size, U_size, variables, variables_length, R_dynamics, M_dynamics, C_dynamics, d_dynamics);
    
    if(status){
        for(i=0;i<U_size;i++){
            for(j=0;j<State_size;j++){
                indexR_i[i*State_size+j] = j;
                indexR_j[i*State_size+j] = i;
            }
        }
    }

    if(status){
        for(i=0;i<S_length-1;i++){
            
            //copy M*S_e_prime for the matrix calculation of m.  These multiplications are performed using BLAS routines.
			cblas_dgemv(CblasRowMajor, CblasNoTrans, m_m, n_m, one, &M_dynamics[(i+1)*State_size*State_size], m_m, &S_e_prime[(i+1)*State_size], p_m, zero, datam, p_m);
            
            //Perform M*S_e_dprime
			cblas_dgemv(CblasRowMajor, CblasNoTrans, m_m, n_m, one, &M_dynamics[(i+1)*State_size*State_size], m_m, &S_e_dprime[(i+1)*State_size], p_m, zero, datac1, p_m);
            
            //C*S_e_prime_s;
			cblas_dgemv(CblasRowMajor, CblasNoTrans, m_m2, n_m2, one, &C_dynamics[(i+1)*State_size*State_size*State_size], n_m2, &S_e_prime_s[(i+1)*State_size*State_size], p_m, zero, datac2, p_m);

            //Now perform all of the assignment;
            //not these indices are all one less than those in matlab;
            if(i==0)
			{
                for(j=0;j<State_size;j++){
                    A_i[j] = j;//c/2
                    A_j[j] = 0;
                    A_v[j] = (datac1[j]+datac2[j]);
					if (A_v[j] * datam[j] >= 0)
					{
						A_i[State_size + j] = j;//m
						A_j[State_size + j] = 1;
					}
					else
					{
						A_i[State_size + j] = j;//m
						A_j[State_size + j] = 1*(2 + U_size) + 1;
					}
                }
                memcpy(&A_v[State_size], datam, State_size*sizeof(double));
                memcpy(&A_i[2*State_size], indexR_i, U_size*State_size*sizeof(int));//-R
                for(j=0;j<U_size*State_size;j++)
				{
                    A_j[2*State_size+j] = indexR_j[j]+2;
                    A_v[2*State_size+j] = R_dynamics[j]*-1;
                }
                A_i[State_size*(2+U_size)]= State_size;//for the -1;
                A_j[State_size*(2+U_size)]= 0;
                A_v[State_size*(2+U_size)]= -1;
                A_i[State_size*(2+U_size)+1]= State_size;//for the 2*dtheta
                A_j[State_size*(2+U_size)+1]= 1;
                A_v[State_size*(2+U_size)+1]= 2*dtheta;
                b_Ax[State_size]=-b_0;
            }
			else
			{
                for(j=0;j<State_size;j++){
                    A_i[i*block_length-State_size+j]= j+i*(State_size+1);//c/2
                    A_j[i*block_length-State_size+j]= i*(2+U_size);                    
					A_v[i*block_length-State_size+j] = (datac1[j] + datac2[j]);
					if ((i == (S_length - 2)) || ((datac1[j] + datac2[j]) * datam[j] >= 0))
					{
						A_i[i*block_length + j] = j + i*(State_size + 1);//m
						A_j[i*block_length + j] = i*(2 + U_size) + 1;
					}
					else
					{
						A_i[i*block_length + j] = j + i*(State_size + 1);//m
						A_j[i*block_length + j] = (i+1)*(2 + U_size) + 1;
					}
                }
                memcpy(&A_v[i*block_length], datam, State_size*sizeof(double));
                for(j=0;j<U_size*State_size;j++){
                    A_i[i*block_length+State_size+j] = indexR_i[j]+i*(State_size+1);
                    A_j[i*block_length+State_size+j] = indexR_j[j]+2+i*(2+U_size);
					A_v[i*block_length + State_size + j] = R_dynamics[i*State_size*U_size+j]*-1;
                }
                A_i[i*block_length+State_size*(1+U_size)]= State_size+i*(State_size+1);//for the -1;
                A_j[i*block_length+State_size*(1+U_size)]= 0+i*(2+U_size);
                A_v[i*block_length+State_size*(1+U_size)]= -1;
                A_i[i*block_length+State_size*(1+U_size)+1]= State_size+i*(State_size+1);//for the 2*dtheta
                A_j[i*block_length+State_size*(1+U_size)+1]= 1+i*(2+U_size);
                A_v[i*block_length+State_size*(1+U_size)+1]= 2*dtheta;
                A_i[i*block_length+State_size*(1+U_size)+2]= State_size+i*(State_size+1);
                A_j[i*block_length+State_size*(1+U_size)+2]=(i-1)*(2+U_size);
                A_v[i*block_length+State_size*(1+U_size)+2]= 1;
                b_Ax[i*(State_size+1)+State_size]=0;
            }
            for(j=0;j<State_size;j++){
				b_Ax[i*(State_size + 1) + j] = d_dynamics[(i+1)*State_size + j] * -1;//the d term
                //if(i==0)
                //    b_Ax[i*(State_size+1)+j]-=b_0*(datac1[j]+datac2[j]);//the c/2*b[0] term.
            }
            
        }
    }
    
    if(status){
        //make sure that the indexed arrays are returned in addition to the entire matrix;
        *pA_i = A_i;
        *pA_j = A_j;
        *pA_v = A_v;
        *pb_Ax = b_Ax;
        
        A.nzmax = block_length*(S_length-1)-State_size;
        A.i = A_i;
        A.p = A_j;
        A.x = A_v;
        A.nz = A.nzmax;
        A.m = (S_length-1)*(State_size+1);
        A.n = (S_length-1)*(2+U_size);
        *pAc = cs_triplet(&A);
        //cs_dupl(Ac);there shouldn't be any duplicates to remove;
        cs_dropzeros(*pAc);
    }
    //free the arrays from the matrix multiply
    if(datam != NULL)
        free(datam);
    if(datac1 != NULL)
        free(datac1);
    if(datac2 != NULL)
        free(datac2);
    
    if(indexR_i != NULL)
        free(indexR_i);
    if(indexR_j != NULL)
        free(indexR_j);

    if(R_dynamics != NULL)
        free(R_dynamics);
    if(M_dynamics != NULL)
        free(M_dynamics);
    if(C_dynamics != NULL)
        free(C_dynamics);
    if(d_dynamics != NULL)
        free(d_dynamics);
    
    if(status)
        return(index_length);
    else
        return -1;
}

//modified from Tim Davis's routine
int so_lu_solve(cs* A, double* b, css** S, int orderOkay) {
    double tol;
    double *x ;
    csn *N ;
    int order, n, ok ;
    
    tol = 1e-10;
    order = 0;
    
    if (!A || !b) return (0) ;		/* check inputs */
    n = A->n ;
    
    if(!orderOkay){
        if(*S!=NULL)
            cs_sfree(*S);
        *S = cs_sqr (A, order, 0);/* ordering and symbolic analysis */
    }
    
    N = cs_lu (A, *S, tol) ;		/* numeric LU factorization */
    x = cs_malloc (n, sizeof (double)) ;
    ok = (*S && N && x) ;
    if (ok)
    {
	cs_ipvec (n, N->Pinv, b, x) ;	/* x = P*b */
	cs_lsolve (N->L, x) ;		/* x = L\x */
	cs_usolve (N->U, x) ;		/* x = U\x */
	cs_ipvec (n, (*S)->Q, x, b) ;	/* b = Q*x */
    }
    cs_free (x) ;
    cs_nfree (N) ;
    if(!ok){
        #ifdef MATLAB_MEX_FILE
                mexErrMsgTxt("LU solve failed");
        #endif
        #ifndef MATLAB_MEX_FILE
                fprintf(stderr, "LU solve failed\n");
        #endif
    }
    return (ok) ;
}

void so_stuffing(double* Abig_v, int* Abig_i, int* Abig_j, double* indexH_v, int* indexH_i, int* indexH_j, int Hindex_length, double* indexA21_v, int* indexA21_i, int* indexA21_j, int A21index_length, int Hc_m, int S_length, int State_size, int U_size){
    int i;
    
    //copy the H array into the big array, but changing ordering to interweave dual varialbes.
    memcpy(Abig_v, indexH_v, Hindex_length*sizeof(double));
    memcpy(Abig_i, indexH_i, Hindex_length*sizeof(int));
    memcpy(Abig_j, indexH_j, Hindex_length*sizeof(int));
    
    
    //copy the A21 array into the big array;
    memcpy(&Abig_v[Hindex_length], indexA21_v, A21index_length*sizeof(double));
    memcpy(&Abig_v[Hindex_length+A21index_length], indexA21_v, A21index_length*sizeof(double));
    memcpy(&Abig_j[Hindex_length], indexA21_j, A21index_length*sizeof(int));
    memcpy(&Abig_i[Hindex_length+A21index_length], indexA21_j, A21index_length*sizeof(int));
    for(i=0;i<A21index_length;i++){
        Abig_i[Hindex_length+i] = indexA21_i[i]+Hc_m;
        Abig_j[Hindex_length+A21index_length+i] = indexA21_i[i]+Hc_m;
    }
}

#ifdef __cplusplus
}
#endif
