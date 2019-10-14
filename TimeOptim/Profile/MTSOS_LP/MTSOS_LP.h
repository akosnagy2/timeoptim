#pragma once

#include "../csparse.h"
#include <stdlib.h>
#include "../MTSOS/MTSOS.h"
/*
typedef struct{
    double kappa;
    double alpha;
    double beta;
    double epsilon;
    int MAX_ITER;
}algorithm_params;

typedef struct{
    double* S;
    int S_length;
    double initial_velocity;
    int State_size;
    int U_size;
}problem_params;

typedef struct{
    int kappa;
    int display;
    int timer;
}algorithm_flags;

typedef struct{
    double* variables;
    int variables_length;
    double* initial_b;
    double* initial_u;
}optional_params;*/

int so_calcHbarrier_LP(double* H_barrier, int S_length, int U_size, int** index_i, int** index_j, double** index_v);
int so_MakeHandG_LP(double* b, double* S_middleArray, double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, int flag_calculate, double** G, cs** pH, int** pindexH_i, int** pindexH_j, double** pindexH_v, double* timers, int istimed, double* variables, int variables_length);
int so_linesearch_infeasible_LP(double* w, double* v, double* x, double* dx, double* b, double* S_middle,double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double beta, double alpha, double* G, cs* A21, double* b_Ax, double* variables, int variables_length);
double so_function_evaluation_LP(double* b, double* S_middle, double* S_prime, double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double* variables, int variables_length);
double so_linesearch_feasible_LP(double* w, double* v, double* x, double* dx, double* b, double* S_middle, double* S_prime,double* S_dprime, double* u, double* a, int S_length, int U_size, int State_size, double dtheta, double kappa, double* G, double alpha, double beta, double fx, double* variables, int variables_length);
int so_MakeA_LP(double* S, double* S_middle, double* S_prime, double* S_dprime, double theta, int U_size, int S_lengt, int State_size, cs** pAc, int** pA_i, int** pA_j, double** pA_v, double** pb_Ax, double b_0, double* variables, int variables_length);
int so_lu_solve_LP(cs* Abig, double* bbig, css** S,int orderOkay);
int so_barrier(const double* S, const double* S_prime, const double* S_dprime, const double* b, const double* a, const double* u, int S_length, int U_size, int State_size, int indicator,double kappa,const double* variables,int variables_length,double* H, double* G, double* F);
int so_dynamics(const double* const S_middle, const double* const S_prime, int S_length, int State_size, int U_size, const double* variables, int variables_length, double* R_dynamics, double* M_dynamics, double* C_dynamics, double* d_dynamics);
void so_stuffing_LP(double* Abig_v, int* Abig_i, int* Abig_j, double* indexH_v, int* indexH_i, int* indexH_j, int Hindex_length, double* indexA21_v, int* indexA21_i, int* indexA21_j, int A21index_length, int Hc_m, int S_length, int State_size, int U_size);
int so_MTSOS_LP(problem_params* p_params, algorithm_flags* flags, algorithm_params* a_params, optional_params* o_params, double** p_b, double** p_u, double** p_v, double* pfx, double** p_timers, int* piterations);