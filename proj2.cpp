#include <iostream>
#include "Matrix.h"
#include "Vector.h"

// Function making a matrix complaint with preoject requirements.
Matrix ourMatrix(unsigned int N, unsigned int e)
{
	Matrix result(N, N);
	double a1 = 5 + e;
	double a2 = -1;
	double a3 = a2;
	for (unsigned int i = 0; i < N; i++)
	{
		result[i][i] = a1;

		if (i + 1 < N)
		{
			result[i][i + 1] = a2;
			result[i + 1][i] = a2;
		}
		if (i + 2 < N)
		{
			result[i][i + 2] = a3;
			result[i + 2][i] = a3;
		}
	}
	return result;
}

// A vector compliant to project requirements.
Vector ourVector(unsigned int N, unsigned int f)
{
	Vector result(N);
	for (unsigned int i = 0; i < N; i++)
		result[i] = sin(i*(f + 1));
	return result;
}


double norm(Vector vector)
{
	double sum = 0;
	for (unsigned int i = 0; i < vector.getLength(); i++)
	{
		sum += vector[i] * vector[i];
	}
	return sqrt(sum);
}


Vector forwardSubst(Matrix L, Vector b)
{
	Vector result(b.getRows());
	unsigned int length = b.getLength();
	for (unsigned int i = 0; i < length; i++)
	{
		result[i] = b[i];
		for (unsigned int j = 0; j < i; j++)
		{
			result[i] -= L[i][j] * result[j];
		}
		result[i] /= L[i][i];
	}
	return result;
}

Vector backwardSubst(Matrix U, Vector b)
{
	Vector result(b.getRows());
	unsigned int length = b.getLength();
	unsigned int cols = U.getCols();
	for (unsigned int i = 0; i < length; i++)
	{
		result[i] = b[i];
		for (unsigned int j = i+1; j < length; j++)
		{
			result[i] -= U[i][j] * result[j];
		}
		result[i] /= U[i][i];
	}
	return result;
}


Matrix LUFactor(Matrix A, Matrix b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	unsigned int m = A.getRows();
	Matrix L(m, m);
	for (unsigned int i = 0; i < m; i++)
		L[i][i] = 1;
	Matrix U = A;
	for (unsigned int k = 0; k < m - 1; k++)
		for (unsigned int j = k + 1; k < m; k++)
		{
			L[j][k] = U[j][k] / U[k][k];
			for (unsigned int l = k; k < m; k++)
				U[j][l] -= (L[j][k] * U[k][l]);
		}

	Matrix test = L*U;
	// For debugging purposes only! (Almost OK, maybe just a precision loss...)
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < m; j++)
		{
			if (test[i][j] != A[i][j])
			{
				std::cout << i << ' ' << j << std::endl;
				std::cout << test[i][j] << ' ' << A[i][j] << std::endl;
			}
		}
	}

	Vector y = forwardSubst(L, b);
	std::cout << norm(L*y - b) << std::endl;
	Vector x = backwardSubst(U, y);
	std::cout << norm(U*x - y) << std::endl;
	return x;
}



Matrix Jacobi(Matrix A, Matrix b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	Matrix L = A;
	Matrix U = A;
	unsigned int rows = A.getRows(), cols = A.getCols();
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j <= i; j++)
		{
			L[i][j] = 0;
			U[rows - i - 1][cols - j - 1] = 0;
		}
	}
	Matrix D = A.diagonal();
	Vector x = b;

	do
	{
		// Implement whatever may be needed! 
	} while (norm(A*x - b) > pow(10, -9));

	// TO BE IMPLEMENTED!
	return Matrix();
}


int main()
{
	Matrix A = ourMatrix(961, 7);
	Vector b = ourVector(961, 5);


	Vector x = LUFactor(A, b);


	system("pause");
	return 0;
}

