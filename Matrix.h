#ifndef MATRIX
#define MATRIX
#include <assert.h>


class Matrix
{
protected:
	unsigned int rows, cols;
	double** arr;

public:

	Matrix();
	Matrix(unsigned int rows, unsigned int cols);
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
	Matrix ones(unsigned int rows, unsigned int cols);
	double& operator() (unsigned int, unsigned int);
	double operator() (unsigned int, unsigned int) const;
	double* operator[] (unsigned int row);
	unsigned int getRows() const;
	unsigned int getCols() const;
};







#endif