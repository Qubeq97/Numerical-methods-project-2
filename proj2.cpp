#include <iostream>
#include "Matrix.h"
#include "Vector.h"




Matrix ourMatrix(unsigned int N, int e)
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

Matrix ourVector(unsigned int N, int f)
{
	Matrix result(N,1);
	for (unsigned int i = 0; i < N; i++)
		result[i][i] = sin(i*(f + 1));
	return result;
}


double norm(Matrix vector)
{
	double sum = 0;
	for (unsigned int i = 0; i < vector.getRows(); i++)
	{
		sum += vector[i][0] * vector[i][0];
	}
	return sqrt(sum);
}

// A must be equal to L * U.
// Not working yet. :(
Vector LUFactor(const Matrix& A, const Vector& b)
{
	unsigned int m = A.getRows();
	Matrix L(m,m);
	for (unsigned int i = 0; i < m; i++)
		L[i][i] = 1;
	Matrix U = A;
	for (unsigned int k=0; k<m-1; k++)
		for (unsigned int j = k + 1; k < m; k++)
		{
			L[j][k] = U[j][k] / U[k][k];
			for (unsigned int l = k-1; k < m; k++)
				U[j][l] -= L[j][j] * U[k][l];
		}
	if (L * U == A)
		std::cout << "YES!\n";
	return Vector(12);
}



// Solves a system of linear equations Ax = b for x using the Jacobi method.
Matrix Jacobi(Matrix A, Matrix b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	// Dividing the A matrix: A = L + U + D
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
	Matrix x = b;
	
	do
	{

	} while (norm(A*x - b) > pow(10, -9));

	return Matrix();
}


int main()
{
	Matrix A = ourMatrix(961, 7);
	Matrix b = ourVector(961, 5);

	Matrix C = A*b;


	Jacobi(A, b);

	return 0;
}

