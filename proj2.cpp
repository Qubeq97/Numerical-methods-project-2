#include <iostream>
#include "Matrix.h"
#include "Vector.h"


// Function making a matrix complaint with preoject requirements.
Matrix ourMatrix(int N, double a1)
{
	Matrix result(N, N);
	double a2 = -1;
	double a3 = a2;
	for (int i = 0; i < N; i++)
	{
		result(i, i) = a1;

		if (i + 1 < N)
		{
			result(i, i + 1) = a2;
			result(i + 1, i) = a2;
		}
		if (i + 2 < N)
		{
			result(i, i + 2) = a3;
			result(i + 2, i) = a3;
		}
	}
	return result;
}

// A vector compliant to project requirements.
Vector ourVector(int N, int f)
{
	Vector result(N);
	for (int i = 0; i < N; i++)
		result(i) = sin(i*(f + 1));
	return result;
}


double norm(const Vector& vector)
{
	double sum = 0;
	for (int i = 0; i < vector.getLength(); i++)
	{
		sum += vector(i) * vector(i);
	}
	return sqrt(sum);
}


Vector forwardSubst(const Matrix& L, const Vector& b)
{
	Vector result(b.getRows());
	int length = b.getLength();
	for (int i = 0; i < length; i++)
	{
		result(i) = b(i);
		for (int j = 0; j < i; j++)
		{
			result(i) -= L(i, j) * result(j);
		}
		result(i) /= L(i, i);
	}
	return result;
}

Vector backwardSubst(const Matrix& U, const Vector& b)
{
	Vector result(b.getRows());
	int length = b.getLength();
	int cols = U.getCols();
	for (int i = length - 1; i >= 0; i--)
	{
		result(i) = b(i);
		for (int j = i + 1; j < length; j++)
		{
			result(i) -= U(i, j) * result(j);
		}
		result(i) /= U(i, i);
	}
	return result;
}


Matrix LUFactor(const Matrix& A, const Matrix& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	int m = A.getRows();
	Matrix L(m, m);
	for (int i = 0; i < m; i++)
		L(i, i) = 1;
	Matrix U = A;
	for (int k = 0; k < m - 1; k++)
		for (int j = k + 1; k < m; k++)
		{
			L(j, k) = U(j, k) / U(k, k);
			for (int l = k; k < m - 1; k++)
				U(j, l) -= (L(j, k) * U(j, l));
		}

	Matrix test = L * U;
	// For debugging purposes. (Almost OK, maybe just a precision loss...)
	for (int i = 0; i < m; i++)
	{
		for (int j = 0; j < m; j++)
		{
			if (test(i, j) != A(i, j))
			{
				std::cout << "LU ERROR: " << i << ' ' << j << std::endl;
				std::cout << test(i, j) << ' ' << A(i, j) << std::endl;
			}
		}
	}

	Vector y = forwardSubst(L, b);
	Vector x = backwardSubst(U, y);
	return x;
}

int Jacobi(const Matrix& A, const Vector& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	Matrix L = A;
	Matrix U = A;
	int rows = A.getRows(), cols = A.getCols();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			L(i, j) = 0;
			U(rows - i - 1, cols - j - 1) = 0;
		}
	}
	Matrix D = A.diagonal();

	Vector x(b.getLength());
	for (int i = 0; i < x.getLength(); i++)
		x(i) = 1;

	L = -L;
	U = -U;


	Matrix sumLU = L + U;

	double resnorm;
	int iters = 0;
	do
	{
		x = forwardSubst(D, sumLU*x) + forwardSubst(D, b);
		resnorm = norm(A*x - b);
		std::cout << resnorm << std::endl;
		iters++;
	} while (resnorm > pow(10, -9));


	return iters;
}

int GaussSeidel(const Matrix& A, const Vector& b)
{
	Matrix L = A;
	Matrix U = A;
	int rows = A.getRows(), cols = A.getCols();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			L(i, j) = 0;
			U(rows - i - 1, cols - j - 1) = 0;
		}
	}
	Matrix D = A.diagonal();

	Vector x(b.getLength());
	for (int i = 0; i < x.getLength(); i++)
		x(i) = 1;

	L = -L;
	U = -U;


	Matrix diffDL = D - L;

	int iters = 0;
	double resnorm;
	do
	{
		x = backwardSubst(diffDL, U*x) + backwardSubst(diffDL, b);
		resnorm = norm(A*x - b);
		std::cout << resnorm << std::endl;
		iters++;
	} while (resnorm > pow(10, -9));

	return iters;
}


int main()
{
	Matrix A = ourMatrix(961, 5 + 7);
	Vector b = ourVector(961, 5);


	int itersa = Jacobi(A, b);
	int itersb = GaussSeidel(A, b);

	std::cout << "Jacobi iters: " << itersa << std::endl;
	std::cout << "Gauss-Seidel iters: " << itersb << std::endl;

	/*A = ourMatrix(961, 3);
	itersa = Jacobi(A, b);
	itersb = GaussSeidel(A, b);
	std::cout << "Jacobi iters: " << itersa << std::endl;
	std::cout << "Gauss-Seidel iters: " << itersb << std::endl;*/

	//In this case, iterative methods DO NOT converge.

	system("pause");
	return 0;
}

