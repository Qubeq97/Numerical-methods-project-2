#include <iostream>
#include <ctime>
#include <fstream>
#include "Matrix.h"
#include "Vector.h"


#define A1 (5+7)
#define F 5

int iters;

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


double LUFactor(const Matrix& A, const Matrix& b)
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
	return norm(A*x - b);
}

double Jacobi(const Matrix& A, const Vector& b)
{
	assert(A.getRows() == A.getCols() && A.getRows() == b.getRows());

	Matrix D = A.diagonal();
	Matrix sumLU = A - D;

	Vector x(b.getLength());
	for (int i = 0; i < x.getLength(); i++)
		x(i) = 1;


	double resnorm;
	iters = 0;
	do
	{
		x = forwardSubst(-D, sumLU*x) + forwardSubst(D, b);
		resnorm = norm(A*x - b);
		iters++;
	} while (resnorm > pow(10, -9));


	return resnorm;
}

double GaussSeidel(const Matrix& A, const Vector& b)
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


	iters = 0;
	double resnorm;
	do
	{
		x = backwardSubst(-sumDL, U*x) + backwardSubst(sumDL, b);
		resnorm = norm(A*x - b);
		iters++;
	} while (resnorm > pow(10, -9));

	return resnorm;
}


int main()
{
	double resnorm;

	Matrix A = ourMatrix(961, A1);
	Vector b = ourVector(961, F);

	clock_t start;

	std::cout << "Test for 961 elements" << std::endl;

	start = clock();
	resnorm = Jacobi(A, b);
	double timeJacobi = double(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Jacobi time: " << timeJacobi << ", iterations: " << iters << ", residuum norm: "<<resnorm<<std::endl;

	start = clock();
	resnorm = GaussSeidel(A, b);
	double timeSeidel = double(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "Gauss-Seidel time: " << timeSeidel << ", iterations: " << iters << ", residuum norm: " << resnorm << std::endl;


	// case C
	A = ourMatrix(961, 3);

	/*itersJ = Jacobi(A, b);
	itersGS = GaussSeidel(A, b);
	std::cout << "Jacobi iters: " << itersJ << std::endl;
	std::cout << "Gauss-Seidel iters: " << itersGS << std::endl;*/
	// Upper ones diverge.

	std::cout << "Test for 961 elements - case C" << std::endl;
	start = clock();
	resnorm = LUFactor(A, b);
	double timeLU = double(clock() - start) / CLOCKS_PER_SEC;
	std::cout << "LU time: " << timeLU << ", residuum norm: " << resnorm << std::endl;


	// In that case, iterative methods DO NOT converge and would iterate infinitely
	// Commented because it would get the program stuck.


	int tab[] = { 100,500,1000,2000,3000,6000 };

	std::ofstream output("results.txt");

	if (!output.is_open())
	{
		std::cerr << "Error opening results file!" << std::endl;
		exit(EXIT_FAILURE);
	}

	output << "N;Jacobi;GS;LU" << std::endl;

	for (int N : tab)
	{
		std::cout << "Calculating for N = " << N << std::endl;

		output << N << ';';

		A = ourMatrix(N, A1);
		b = ourVector(N, F);


		start = clock();
		resnorm = Jacobi(A, b);
		timeJacobi = double(clock() - start) / CLOCKS_PER_SEC;
		output << timeJacobi << ';';

		start = clock();
		resnorm = GaussSeidel(A, b);
		timeSeidel = double(clock() - start) / CLOCKS_PER_SEC;
		output << timeSeidel << ';';

		start = clock();
		resnorm = LUFactor(A, b);
		double timeLU = double(clock() - start) / CLOCKS_PER_SEC;
		output << timeLU<< std::endl;

	}

	output.close();

	std::cout << "Results calculated successfully." << std::endl;
	system("pause");
	return 0;
}

