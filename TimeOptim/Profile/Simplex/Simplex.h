#pragma once

#define SIMPLEX_M (100)
#define SIMPLEX_N (100)

static const double epsilon = 1.0e-8;

typedef struct 
{
	int m, n; // m=rows, n=columns, mat[m x n]
	double mat[SIMPLEX_M][SIMPLEX_N];
} Tableau;

void simplex(Tableau *tab);
void get_optimal_vector(Tableau *tab, double* b);