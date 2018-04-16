#ifndef MATRIX
#define MATRIX
#include <assert.h>


class Matrix
{
protected:
	int rows, cols;
	double** arr;

public:

	Matrix();
	Matrix(int rows, int cols);
	~Matrix();
	Matrix(const Matrix& m);
	Matrix operator=(const Matrix& m);
	Matrix operator=(Matrix&& m);
	Matrix operator+(const Matrix& m) const;
	Matrix operator-(const Matrix& m) const;
	Matrix operator*(const Matrix& m) const;
	Matrix operator-() const;
	bool operator == (const Matrix& m) const;
	Matrix diagonal() const;
	Matrix transposed();
	Matrix ones(int rows, int cols);
	double& operator() (int, int);
	double operator() (int, int) const;
	double* operator[] (int row);
	int getRows() const;
	int getCols() const;
};







#endif