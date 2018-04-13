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
	Matrix operator+(const Matrix& m);
	Matrix operator-(const Matrix& m);
	Matrix operator*(const Matrix& m);
	bool operator == (const Matrix& m);
	Matrix transposed();
	Matrix ones(unsigned int rows, unsigned int cols);
	double* operator[] (unsigned int row);
	unsigned int getRows() const;
	unsigned int getCols() const;
};

Matrix diagonal(Matrix & m);





#endif