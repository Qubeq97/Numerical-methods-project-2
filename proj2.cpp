#include <iostream>
#include <ctime>
#include "Matrix.h"
#include "Vector.h"

#define A1 (5+7)
#define F 5

// Function making a matrix complaint with project requirements.
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


Vector LUFactor(const Matrix& A, const Matrix& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());
	int m = A.getRows();
	Matrix L(m, m);
	for (int i = 0; i < m; i++)
		L(i, i) = 1;
	Matrix U = A;
	for (int k = 0; k < m - 1; k++)
		for (int j = k + 1; j < m; j++)
		{
			L(j, k) = U(j, k) / U(k, k);
			for (int l = k; l < m; l++)
				U(j, l) -= (L(j, k) * U(k, l));
		}


	Vector y = forwardSubst(L, b);
	Vector x = backwardSubst(U, y);
	return x;
}

int Jacobi(const Matrix& A, const Vector& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());

	Matrix D = A.diagonal();
	Matrix sumLU = A - D;

	Vector x(b.getLength());
	for (int i = 0; i < x.getLength(); i++)
		x(i) = 1;


	double resnorm;
	int iters = 0;
	do
	{
		x = forwardSubst(-D, sumLU*x) + forwardSubst(D, b);
		resnorm = norm(A*x - b);
		iters++;
	} while (resnorm > pow(10, -9));


	return iters;
}

int GaussSeidel(const Matrix& A, const Vector& b)
{
	Matrix U = A;
	int rows = A.getRows(), cols = A.getCols();
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j <= i; j++)
		{
			U(rows - i - 1, cols - j - 1) = 0;
		}
	}

	Matrix sumDL = A - U;

	Vector x(b.getLength());
	for (int i = 0; i < x.getLength(); i++)
		x(i) = 1;


	int iters = 0;
	double resnorm;
	do
	{
		x = backwardSubst(-sumDL, U*x) + backwardSubst(sumDL, b);
		resnorm = norm(A*x - b);
		iters++;
	} while (resnorm > pow(10, -9));

	return iters;
}


int main()
{
	Matrix A = ourMatrix(961, A1);
	Vector b = ourVector(961, F);

	clock_t start;

	std::cout << "Test for 961 elements" << std::endl;

	start = clock();
	int itersJ = Jacobi(A, b);
	double timeJacobi = double(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Jacobi time: " << timeJacobi << ", iterations: " << itersJ << std::endl;

	start = clock();
	int itersGS = GaussSeidel(A, b);
	double timeSeidel = double(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Gauss-Seidel time: " << timeSeidel << ", iterations: " << itersGS << std::endl;


	Vector luResult = LUFactor(A, b);
	std::cout << "Norm of LU result residuum: " << norm(A*luResult - b) << std::endl;

	/*A = ourMatrix(961, 3);
	itersa = Jacobi(A, b);
	itersb = GaussSeidel(A, b);
	std::cout << "Jacobi iters: " << itersa << std::endl;
	std::cout << "Gauss-Seidel iters: " << itersb << std::endl;*/

	// In that case, iterative methods DO NOT converge and would iterate infinitely
	// Commented because it would get the program stuck.


	int tab[] = { 100,500,1000,2000,3000,6000 };

	for (int N : tab)
	{
		std::cout << "Test for " << N << " elements" << std::endl;
		A = ourMatrix(N, A1);
		b = ourVector(N, F);


		start = clock();
		itersJ = Jacobi(A, b);
		timeJacobi = double(clock() - start) / CLOCKS_PER_SEC;
		std::cout << "Jacobi time: " << timeJacobi << ", iterations: " << itersJ << std::endl;

		start = clock();
		itersGS = GaussSeidel(A, b);
		timeSeidel = double(clock() - start) / CLOCKS_PER_SEC;
		std::cout << "Gauss-Seidel time: " << timeSeidel << ", iterations: " << itersGS << std::endl;

	}

	system("pause");
	return 0;
}

