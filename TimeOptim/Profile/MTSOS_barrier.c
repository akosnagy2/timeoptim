#include <math.h>
//#include "MTSOS.h"
#ifdef MATLAB_MEX_FILE
#include "blas.h" //to ensure that the blas library is included append the argument -lmwblas to the compile
#endif
#ifndef MATLAB_MEX_FILE
#include <cblas.h>
#endif
/*A barrier for the front wheel drive friction circle car*/


int so_barrier(const double* S, const double* S_prime, const double* S_dprime, const double* b, const double* a, const double* u, int S_length, int U_size, int State_size, int indicator, double kappa, const double* variables, int variables_length, double* H_barrier, double* G_barrier, double* F_barrier){
	//F, G, and H are all passed in as NULL pointers;
	double *fxA1, *fxA2, *fxV, *fxD1Max, *fxD1Min, *fxD2Max, *fxD2Min, *fxD3Max, *fxD3Min;
	int valid, block_length;
	int i, status, var_offset;

	status = 1;//a variable to make sure memory allocation runs smoothly

	//calculate the value of f(x)-c for a constraint f(x) < c;
	fxA1 = malloc((S_length - 1)*sizeof(double));
	fxA2 = malloc((S_length - 1)*sizeof(double));
	fxV = malloc((S_length - 1)*sizeof(double));
	fxD1Max = malloc((S_length - 1)*sizeof(double));
	fxD1Min = malloc((S_length - 1)*sizeof(double));
	fxD2Max = malloc((S_length - 1)*sizeof(double));
	fxD2Min = malloc((S_length - 1)*sizeof(double));
	fxD3Max = malloc((S_length - 1)*sizeof(double));
	fxD3Min = malloc((S_length - 1)*sizeof(double));

	if (fxA1 == NULL || fxA2 == NULL || fxV == NULL || fxD1Max == NULL || fxD1Min == NULL || fxD2Max == NULL || fxD2Min == NULL || fxD3Max == NULL || fxD3Min == NULL)
	{//insufficient memory
		status = 0;
	}
	else
	{
		valid = 1;//a parameter to ensure that the inequality constraint is always met.
		for (i = 0; i<(S_length - 1) && valid == 1; i++)
		{
			var_offset = variables_length - S_length - (S_length - 1);

			/* Acceleration constraints */
			fxA1[i] = a[i] - variables[var_offset + i];
			fxA2[i] = -a[i] - variables[var_offset + i];
			
			/* Velocity constraints */
			fxV[i] = b[i] - variables[var_offset + S_length + i];

			/* Dynamics constrains */
			if (i == (S_length - 2))
			{
				fxD1Max[i] = -10;
				fxD2Max[i] = -10;
				fxD3Max[i] = -10;
				fxD1Min[i] = -10;
				fxD2Min[i] = -10;
				fxD3Min[i] = -10;
			}
			else
			{
				fxD1Max[i] = u[i*U_size] - variables[var_offset - 3];
				fxD2Max[i] = u[i*U_size + 1] - variables[var_offset - 2];
				fxD3Max[i] = u[i*U_size + 2] - variables[var_offset - 1];
				fxD1Min[i] = -u[i*U_size] - variables[var_offset - 3];
				fxD2Min[i] = -u[i*U_size + 1] - variables[var_offset - 2];
				fxD3Min[i] = -u[i*U_size + 2] - variables[var_offset - 1];
			}

			if (fxA1[i]>0 || fxA2[i]>0 || fxV[i]>0 || fxD1Max[i]>0 || fxD2Max[i]>0 || fxD3Max[i]>0 || fxD1Min[i]>0 || fxD2Min[i]>0 || fxD3Min[i]>0)
				valid = 0;
		}		

		//if the constraints are violated, then the NULL pointers are returned;
		if (valid != 0)
		{
			if (indicator < 3)
			{
				if (indicator < 2 && status)
				{
					//Allocate memory to store an array of Hessians(the first (2+U_size)*(2+U_size) block refers to the first Hessian)
					block_length = (2 + U_size)*(2 + U_size);//the size of one Hessian block
					//zero the allocated memory (for speed up this could be incorporated into the allocation loop as it is with G.
					for (i = 0; i<(2 + U_size)*(2 + U_size)*(S_length - 1); i++)
						H_barrier[i] = 0;
					for (i = 0; i<(S_length - 1) && status; i++)
					{
						H_barrier[i*block_length + (2 + U_size) * 0 + 0] = kappa*(1 / (fxV[i] * fxV[i]));
						H_barrier[i*block_length + (2 + U_size) * 1 + 1] = kappa*(1 / (fxA1[i] * fxA1[i]) + 1 / (fxA2[i] * fxA2[i]));
						H_barrier[i*block_length + (2 + U_size) * 2 + 2] = kappa*(1 / (fxD1Max[i] * fxD1Max[i]) + 1 / (fxD1Min[i] * fxD1Min[i]));
						H_barrier[i*block_length + (2 + U_size) * 3 + 3] = kappa*(1 / (fxD2Max[i] * fxD2Max[i]) + 1 / (fxD2Min[i] * fxD2Min[i]));
						H_barrier[i*block_length + (2 + U_size) * 4 + 4] = kappa*(1 / (fxD3Max[i] * fxD3Max[i]) + 1 / (fxD3Min[i] * fxD3Min[i]));
					}
				}
				if ((indicator == 0 || indicator == 2) && status)
				{
					for (i = 0; i<S_length - 1; i++)
					{
						G_barrier[(2 + U_size)*i] = -kappa*(1 / fxV[i]);
						G_barrier[(2 + U_size)*i + 1] = -kappa*(1 / fxA1[i] - 1 / fxA2[i]);
						G_barrier[(2 + U_size)*i + 2] = -kappa*(1 / fxD1Max[i] - 1 / fxD1Min[i]);
						G_barrier[(2 + U_size)*i + 3] = -kappa*(1 / fxD2Max[i] - 1 / fxD2Min[i]);
						G_barrier[(2 + U_size)*i + 4] = -kappa*(1 / fxD3Max[i] - 1 / fxD3Min[i]);
					}
				}
			}
			else
			{
				for (i = 0; i<S_length - 1; i++)
					F_barrier[i] = -kappa*(log(-fxA1[i]) + log(-fxA2[i]) + log(-fxV[i]) + log(-fxD1Max[i]) + log(-fxD2Max[i]) + log(-fxD3Max[i]) + log(-fxD1Min[i]) + log(-fxD2Min[i]) + log(-fxD3Min[i]));
			}
		}
	}
	if (fxA1 != NULL)
		free(fxA1);
	if (fxA2 != NULL)
		free(fxA2);
	if (fxV != NULL)
		free(fxV);
	if (fxD1Max != NULL)
		free(fxD1Max);
	if (fxD2Max != NULL)
		free(fxD2Max);
	if (fxD3Max != NULL)
		free(fxD3Max);
	if (fxD1Min != NULL)
		free(fxD1Min);
	if (fxD2Min != NULL)
		free(fxD2Min);
	if (fxD3Min != NULL)
		free(fxD3Min);
	if (status)
		return valid;
	else
		return -1;
}
