#ifndef MATRIX
#define MATRIX
#include <assert.h>


class Matrix
{
protected:
	int rows, cols;
	double** arr=nullptr;

public:

	Matrix();
	Matrix(int rows, int cols);
	~Matrix();
	Matrix(const Matrix& m);
	Matrix(Matrix&& m);
	Matrix operator=(const Matrix& m);
	Matrix operator=(Matrix&& m);
	Matrix operator+(const Matrix& m) const;
	Matrix operator-(const Matrix& m) const;
	Matrix operator*(const Matrix& m) const;
	Matrix operator-() const;
	bool operator == (const Matrix& m) const;
	double& operator() (int, int);
	double operator() (int, int) const;
	int getRows() const;
	int getCols() const;
};

#endif