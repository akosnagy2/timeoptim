//#include "MTSOS.h"
#include <math.h>

int so_dynamics(const double* const S, const double* const S_prime, int S_length, int State_size, int U_size, const double* variables, int variables_length, double* R_dynamics, double* M_dynamics, double* C_dynamics, double* d_dynamics){
	int i,j,k,block_length;

	double Ix1, Ix2, Ix3, Iy1, Iy2, Iy3, Iz1, Iz2, Iz3;
	double m1, m2, m3;
	double r0, r1, r2;
	double l0, l1, l2;
	double C[3][3][3];

	for (i = 0; i < 3; ++i)
		for (j = 0; j < 3; ++j)
			for (k = 0; k < 3; ++k)
				C[i][j][k] = 0;

	//Ix
	Ix1 = variables[0];
	Ix2 = variables[1];
	Ix3 = variables[2];
	//Iy
	Iy1 = variables[3];
	Iy2 = variables[4];
	Iy3 = variables[5];
	//Iz
	Iz1 = variables[6];
	Iz2 = variables[7];
	Iz3 = variables[8];
	//m
	m1 = variables[9];
	m2 = variables[10];
	m3 = variables[11];
	//r
	r0 = variables[12];
	r1 = variables[13];
	r2 = variables[14];
	//l
	l0 = variables[15];
	l1 = variables[16];
	l2 = variables[17];

	block_length = U_size*State_size;
	for (i = 0; i<S_length; i++)
		for (j = 0; j < block_length; ++j)
			if (j % 4 ==  0)
				R_dynamics[i*block_length + j] = 1;
			else
				R_dynamics[i*block_length + j] = 0;	

	block_length = State_size*State_size;
	for (i = 0; i < S_length; i++)
	{
		double s2 = sin(S[i*State_size + 1]);
		double c2 = cos(S[i*State_size + 1]);
		double c3 = cos(S[i*State_size + 2]);
		double s23 = sin(S[i*State_size + 1] + S[i*State_size + 2]);
		double c23 = cos(S[i*State_size + 1] + S[i*State_size + 2]);

		//First row
		M_dynamics[i*block_length + 0] = Iy2*s2*s2 + Iy3*s23*s23 + Iz1 + Iz2*c2*c2 + Iz3*c23*c23 + m2*r1*r1*c2*c2 + m3*(l1*c2 + r2*c23)*(l1*c2 + r2*c23);
		M_dynamics[i*block_length + 1] = 0;
		M_dynamics[i*block_length + 2] = 0;
		//Second row
		M_dynamics[i*block_length + 3] = 0;
		M_dynamics[i*block_length + 4] = Ix2 + Ix3 + m3*l1*l1 + m2*r1*r1 + m3*r2*r2 + 2 * m3*l1*r2*c3;
		M_dynamics[i*block_length + 5] = Ix3 + m3*r2*r2 + m3*l1*r2*c3;
		//Third row
		M_dynamics[i*block_length + 6] = 0;
		M_dynamics[i*block_length + 7] = Ix3 + m3*r2*r2 + m3*l1*r2*c3;
		M_dynamics[i*block_length + 8] = Ix3 + m3*r2*r2;
	}

	block_length = State_size*State_size*State_size;
	for (i = 0; i < S_length; i++)
	{
		double s2 = sin(S[i*State_size + 1]);
		double s3 = sin(S[i*State_size + 2]);
		double c2 = cos(S[i*State_size + 1]);
		double c3 = cos(S[i*State_size + 2]);
		double s23 = sin(S[i*State_size + 1] + S[i*State_size + 2]);
		double c23 = cos(S[i*State_size + 1] + S[i*State_size + 2]);

		//111-113
		C[0][0][1] = (Iy2 - Iz2 - m2*r1*r1)*c2*s2 + (Iy3 - Iz3)*c23*s23 - m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
		C[0][0][2] = (Iy3 - Iz3)*c23*s23 - m3*r2*s23*(l1*c2 + r2*c23);
		//121-123 
		C[0][1][0] = (Iy2 - Iz2 - m2*r1*r1)*c2*s2 + (Iy3 - Iz3)*c23*s23 - m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
		//131-133
		C[0][2][0] = (Iy3 - Iz3)*c23*s23 - m3*r2*s23*(l1*c2 + r2*c23);
		//211-213
		C[1][0][0] = (Iz2 - Iy2 + m2*r1*r1)*c2*s2 + (Iz3 - Iy3)*c23*s23 + m3*(l1*c2 + r2*c23)*(l1*s2 + r2*s23);
		//221-223
		C[1][1][2] = -l1*m3*r2*s3;
		//231-233
		C[1][2][1] = -l1*m3*r2*s3;
		C[1][2][2] = -l1*m3*r2*s3;
		//311-313
		C[2][0][0] = (Iz3 - Iy3)*c23*s23 + m3*r2*s23*(l1*c2 + r2*c23);
		//321-323
		C[2][1][1] = l1*m3*r2*s3;
	
		//First row;
		C_dynamics[i*block_length + 0] = C[0][0][0];
		C_dynamics[i*block_length + 1] = C[0][1][0];
		C_dynamics[i*block_length + 2] = C[0][2][0];
		C_dynamics[i*block_length + 3] = C[0][0][1];
		C_dynamics[i*block_length + 4] = C[0][1][1];
		C_dynamics[i*block_length + 5] = C[0][2][1];
		C_dynamics[i*block_length + 6] = C[0][0][2];
		C_dynamics[i*block_length + 7] = C[0][1][2];
		C_dynamics[i*block_length + 8] = C[0][2][2];
		//Second row;
		C_dynamics[i*block_length + 9] = C[1][0][0];
		C_dynamics[i*block_length + 10] = C[1][1][0];
		C_dynamics[i*block_length + 11] = C[1][2][0];
		C_dynamics[i*block_length + 12] = C[1][0][1];
		C_dynamics[i*block_length + 13] = C[1][1][1];
		C_dynamics[i*block_length + 14] = C[1][2][1];
		C_dynamics[i*block_length + 15] = C[1][0][2];
		C_dynamics[i*block_length + 16] = C[1][1][2];
		C_dynamics[i*block_length + 17] = C[1][2][2];
		//Third row;
		C_dynamics[i*block_length + 18] = C[2][0][0];
		C_dynamics[i*block_length + 19] = C[2][1][0];
		C_dynamics[i*block_length + 20] = C[2][2][0];
		C_dynamics[i*block_length + 21] = C[2][0][1];
		C_dynamics[i*block_length + 22] = C[2][1][1];
		C_dynamics[i*block_length + 23] = C[2][2][1];
		C_dynamics[i*block_length + 24] = C[2][0][2];
		C_dynamics[i*block_length + 25] = C[2][1][2];
		C_dynamics[i*block_length + 26] = C[2][2][2];
	}

	block_length = State_size;
	for (i = 0; i < S_length; i++)
	{
		double c2 = cos(S[i*State_size + 1]);
		double c23 = cos(S[i*State_size + 1] + S[i*State_size + 2]);
		d_dynamics[i*block_length + 0] = 0;
		d_dynamics[i*block_length + 1] = -(m2*9.81*r1 + m3*9.81*l1)*c2 - m3*9.81*r2*c23;
		d_dynamics[i*block_length + 2] = -m3*9.81*r2*c23;
	}

	return 1;
}
