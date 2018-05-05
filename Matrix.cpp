#include "Matrix.h"
#include <stdlib.h>
#include <algorithm>

Matrix::Matrix()
{
	rows = 0;
	cols = 0;
	arr = nullptr;
}

Matrix::Matrix(int rows, int cols)
{
	this->rows = rows;
	this->cols = cols;
	arr = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (int j = 0; j < cols; j++)
			arr[i][j] = 0;
	}
}


Matrix::~Matrix()
{
	if (arr != nullptr)
	{
		for (int i = 0; i < rows; i++)
			delete[] arr[i];
		delete[] arr;
	}
	arr = nullptr;
}

Matrix::Matrix(const Matrix& m)
{
	rows = m.rows;
	cols = m.cols;
	arr = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (int j = 0; j < cols; j++)
			arr[i][j] = m.arr[i][j];
	}
}

Matrix::Matrix(Matrix && m)
{
	std::swap(rows, m.rows);
	std::swap(cols, m.cols);
	std::swap(arr, m.arr);
}



Matrix Matrix:: operator=(const Matrix& m)
{
	rows = m.rows;
	cols = m.cols;
	arr = new double*[rows];
	for (int i = 0; i < rows; i++)
	{
		arr[i] = new double[cols];
		for (int j = 0; j < cols; j++)
			arr[i][j] = m.arr[i][j];
	}
	return *this;
}

Matrix Matrix:: operator=(Matrix&& m)
{
	std::swap(rows, m.rows);
	std::swap(cols, m.cols);
	std::swap(arr, m.arr);
	return *this;
}

Matrix Matrix:: operator+(const Matrix& m) const
{
	Matrix result(*this);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			result.arr[i][j] += m.arr[i][j];
	}
	return result;
}

Matrix Matrix:: operator-(const Matrix& m) const
{
	Matrix result(*this);
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			result.arr[i][j] -= m.arr[i][j];
	}
	return result;
}

Matrix Matrix::operator*(const Matrix& m) const
{
	Matrix result(this->rows, m.cols);
	for (int i = 0; i < this->rows; i++)
	{
		for (int j = 0; j < m.cols; j++)
		{
			double sum = 0;
			for (int k = 0; k < this->cols; k++)
			{
				sum += this->arr[i][k] * m.arr[k][j];
			}
			result.arr[i][j] = sum;
		}
	}
	return result;
}

Matrix Matrix::operator-() const
{
	Matrix result(*this);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			result.arr[i][j] *= (-1);
	return result;
}

bool Matrix::operator == (const Matrix& m) const
{
	if (arr == nullptr || m.arr == nullptr)
		return false;
	if (cols != m.cols || rows != m.rows)
		return false;
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++)
			if (arr[i][j] != m.arr[i][j])
				return false;
	return true;
}



double& Matrix::operator()(int row, int col)
{
	assert(row < rows && col < cols);
	return arr[row][col];
}

double Matrix::operator()(int row, int col) const
{
	assert(row < rows && col < cols);
	return arr[row][col];
}



int Matrix::getRows() const
{
	return rows;
}
int Matrix::getCols() const
{
	return cols;
}
