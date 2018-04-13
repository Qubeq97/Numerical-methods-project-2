#include "Matrix.h"
#include <stdlib.h>
#include <algorithm>

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	arr = nullptr;
}

Matrix::Matrix(unsigned int rows, unsigned int cols)
	:rows(rows), cols(cols)
{
	arr = (double**)calloc(rows, sizeof(double*));
	for (unsigned int i = 0; i < rows; i++)
	{
		arr[i] = (double*)calloc(cols, sizeof(double));
	}
}


Matrix::~Matrix()
{
	if (arr != nullptr)
	{
		for (unsigned int i = 0; i < rows; i++)
			free(arr[i]);
		free(arr);
	}
}

Matrix::Matrix(const Matrix& m)
{
	rows = m.rows;
	cols = m.cols;
	arr = new double*[rows];
	for (unsigned int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (unsigned int j = 0; j < cols; j++)
			arr[i][j] = m.arr[i][j];
	}
}

Matrix Matrix:: operator=(const Matrix& m)
{
	rows = m.rows;
	cols = m.cols;
	arr = new double*[rows];
	for (unsigned int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (unsigned int j = 0; j < cols; j++)
			arr[i][j] = m.arr[i][j];
	}
	return *this;
}

Matrix Matrix:: operator=(Matrix&& m)
{
	rows = m.rows;
	cols = m.cols;
	std::swap(arr, m.arr);
	return *this;
}

Matrix Matrix:: operator+(const Matrix& m)
{
	Matrix result(*this);
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
			result.arr[i][j] += m.arr[i][j];
	}
	return result;
}

Matrix Matrix:: operator-(const Matrix& m)
{
	Matrix result(*this);
	for (unsigned int i = 0; i < rows; i++)
	{
		for (unsigned int j = 0; j < cols; j++)
			result.arr[i][j] -= m.arr[i][j];
	}
	return result;
}

Matrix Matrix::operator*(const Matrix& m)
{
	Matrix result(this->rows, m.cols);
	for (unsigned int i = 0; i < this->rows; i++)
	{
		for (unsigned int j = 0; j < m.cols; j++)
		{
			double sum = 0;
			for (unsigned int k = 0; k < this->cols; k++)
			{
				sum += this->arr[i][k] * m.arr[k][j];
			}
			result.arr[i][j] = sum;
		}
	}
	return result;
}

bool Matrix:: operator == (const Matrix& m)
{
	if (arr == nullptr || m.arr == nullptr)
		return false;
	if (cols != m.cols || rows != m.rows)
		return false;
	for (unsigned int i = 0; i < rows; i++)
		for (unsigned int j = 0; j < cols; j++)
			if (arr[i][j] != m.arr[i][j])
				return false;
	return true;
}

Matrix Matrix::transposed()
{
	Matrix result(cols, rows);
	for (unsigned int i = 0; i < cols; i++)
		for (unsigned int j = 0; j < rows; j++)
			result.arr[i][j] = this->arr[j][i];
	return result;
}

Matrix Matrix::ones(unsigned int rows, unsigned int cols)
{
	Matrix result;
	result.arr = arr = (double**)calloc(rows, sizeof(double*));
	for (unsigned int i = 0; i < rows; i++)
	{
		arr[i] = (double*)calloc(cols, sizeof(double));
		arr[i][i] = 1;
	}
	return result;
}



double* Matrix::operator[] (unsigned int row)
{
	return arr[row];
}


unsigned int Matrix::getRows() const
{
	return rows;
}
unsigned int Matrix::getCols() const
{
	return cols;
}



Matrix diagonal(Matrix & m)
{
	assert(m.getCols() == m.getRows());
	Matrix result = Matrix(m.getCols(), m.getRows());
	for (unsigned int i = 0; i < m.getRows(); i++)
		result[i][i] = m[i][i];
	return result;
}