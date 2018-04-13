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

Vector ourVector(unsigned int N, int f)
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
			for (unsigned int l = k; k < m; k++)
				U[j][l] -= L[j][j] * U[k][l];
		}
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < m; j++)
		{
			if (L[i][j] != 0 || U[i][j] != 0)
			{
				std::cout << "NONZERO!" << std::endl;
				return Vector(12);
			}
		}
	}
	return Vector(12);
}

int main()
{
	Matrix A = ourMatrix(961, 7);
	Vector b = ourVector(961, 5);


	Matrix C = A*b;


	LUFactor(A, b);

	return 0;
}

