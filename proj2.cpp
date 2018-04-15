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
		result(i,i) = a1;

		if (i + 1 < N)
		{
			result(i,i + 1) = a2;
			result(i + 1,i) = a2;
		}
		if (i + 2 < N)
		{
			result(i,i + 2) = a3;
			result(i + 2,i) = a3;
		}
	}
	return result;
}

// A vector compliant to project requirements.
Vector ourVector(unsigned int N, unsigned int f)
{
	Vector result(N);
	for (unsigned int i = 0; i < N; i++)
		result(i) = sin(i*(f + 1));
	return result;
}


double norm(const Vector& vector)
{
	double sum = 0;
	for (unsigned int i = 0; i < vector.getLength(); i++)
	{
		sum += vector(i) * vector(i);
	}
	return sqrt(sum);
}


Vector forwardSubst(const Matrix& L, const Vector& b)
{
	Vector result(b.getRows());
	unsigned int length = b.getLength();
	for (unsigned int i = 0; i < length; i++)
	{
		result(i) = b(i);
		for (unsigned int j = 0; j < i; j++)
		{
			result(i) -= L(i,j) * result(j);
		}
		result(i) /= L(i,i);
	}
	return result;
}

Vector backwardSubst(const Matrix& U, const Vector& b)
{
	Vector result(b.getRows());
	unsigned int length = b.getLength();
	unsigned int cols = U.getCols();
	for (unsigned int i = 0; i < length; i++)
	{
		result(i) = b(i);
		for (unsigned int j = i + 1; j < length; j++)
		{
			result(i) -= U(i,j) * result(j);
		}
		result(i) /= U(i,i);
	}
	return result;
}


Matrix LUFactor(const Matrix& A, const Matrix& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	unsigned int m = A.getRows();
	Matrix L(m, m);
	for (unsigned int i = 0; i < m; i++)
		L(i,i) = 1;
	Matrix U = A;
	for (unsigned int k = 0; k < m - 1; k++)
		for (unsigned int j = k + 1; k < m; k++)
		{
			L(j,k) = U(j,k) / U(k,k);
			for (unsigned int l = k; k < m - 1; k++)
				U(j,l) -= (L(j,k) * U(j,l));
		}


	Matrix test = L*U;
	// For debugging purposes only! (Almost OK, maybe just a precision loss...)
	for (unsigned int i = 0; i < m; i++)
	{
		for (unsigned int j = 0; j < m; j++)
		{
			if (test(i,j) != A(i,j))
			{
				std::cout << "LU ERROR: " << i << ' ' << j << std::endl;
				std::cout << test(i,j) << ' ' << A(i,j) << std::endl;
			}
		}
	}

	Vector y = forwardSubst(L, b);
	Vector x = backwardSubst(U, y);
	return x;
}


Vector Jacobi(const Matrix& A, const Vector& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	Matrix L = A;
	Matrix U = A;
	unsigned int rows = A.getRows(), cols = A.getCols();
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j <= i; j++)
		{
			L(i,j) = 0;
			U(rows - i - 1,cols - j - 1) = 0;
		}
	}
	Matrix D = A.diagonal();

	Vector x(b.getLength());
	for (unsigned int i = 0; i < x.getLength(); i++)
		x(i) = 1;

	L = -L;
	U = -U;
	assert(A == -L -U +D);

	Matrix sumLU = L + U;

	double resnorm;
	do
	{
		x = forwardSubst(D, sumLU*x) + forwardSubst(D, b);
		resnorm = norm(A*x - b);
		std::cout << resnorm << std::endl;
	} while (resnorm > pow(10, -9));


	return x;
}

Vector GaussSeidel(const Matrix& A, const Vector& b)
{
	Matrix L = A;
	Matrix U = A;
	unsigned int rows = A.getRows(), cols = A.getCols();
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j <= i; j++)
		{
			L(i,j) = 0;
			U(rows-i-1,cols-j-1) = 0;
		}
	}
	Matrix D = A.diagonal();

	Vector x(b.getLength());
	for (unsigned int i = 0; i < x.getLength(); i++)
		x(i) = 1;

	L = -L;
	U = -U;
	assert(A == -L -U +D);

	Matrix diffDL = D - L;

	double resnorm;
	do
	{
		x = forwardSubst(diffDL, U*x) + forwardSubst(diffDL, b);
		resnorm = norm(A*x - b);
		std::cout << resnorm << std::endl;
	} while (resnorm > pow(10, -9));
	return x;
}


int main()
{
#ifndef OMGWTF
	Matrix A = ourMatrix(961, 7);
	Vector b = ourVector(961, 5);
#else
	Matrix A(3, 3);
	A[0][0] = 2;
	A[0][1] = -2;
	A[0][2] = -2;
	A[1][0] = 5;
	A[1][1] = 2;
	A[1][2] = 3;
	A[2][0] = -1;
	A[2][1] = 3;
	A[2][2] = 4;

	Vector b(3);
	b[0] = -2;
	b[1] = 8;
	b[2] = 4;
#endif

	Vector x = LUFactor(A, b);
	std::cout << norm(A*x - b) << std::endl;


	system("pause");
	return 0;
}

