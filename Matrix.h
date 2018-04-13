#ifndef MATRIX
#define MATRIX
#include <assert.h>

class Matrix
{
protected:
	unsigned int rows, cols;
	double** arr;

public:

	Matrix()
	{
		rows = 0;
		cols = 0;
		arr = nullptr;
	}

	Matrix(unsigned int rows, unsigned int cols)
		:rows(rows), cols(cols)
	{
		arr = (double**)calloc(rows, sizeof(double*));
		for (unsigned int i = 0; i < rows; i++)
		{
			arr[i] =(double*) calloc(cols, sizeof(double));
		}
	}


	~Matrix()
	{
		if (arr != nullptr)
		{
			for (unsigned int i = 0; i < rows; i++)
				free(arr[i]);
			free(arr);
		}
	}

	Matrix(const Matrix& m)
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

	Matrix operator=(const Matrix& m)
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

	Matrix operator=(Matrix&& m)
	{
		rows = m.rows;
		cols = m.cols;
		std::swap(arr, m.arr);
		return *this;
	}

	Matrix operator+(const Matrix& m)
	{
		Matrix result(*this);
		for (unsigned int i = 0; i < rows; i++)
		{
			for (unsigned int j = 0; j < cols; j++)
				result.arr[i][j] += m.arr[i][j];
		}
		return result;
	}

	Matrix operator-(const Matrix& m)
	{
		Matrix result(*this);
		for (unsigned int i = 0; i < rows; i++)
		{
			for (unsigned int j = 0; j < cols; j++)
				result.arr[i][j] -= m.arr[i][j];
		}
		return result;
	}

	Matrix operator*(const Matrix& m)
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

	Matrix transposed(const Matrix& m)
	{
		Matrix result(cols, rows);
		for (unsigned int i = 0; i < cols; i++)
			for (unsigned int j = 0; j < rows; j++)
				result.arr[i][j] = this->arr[j][i];
		return result;
	}

	Matrix ones(unsigned int rows, unsigned int cols)
	{
		Matrix result;
		result.arr = arr = (double**)calloc(rows, sizeof(double*));
		for (unsigned int i = 0; i < rows; i++)
		{
			arr[i] = (double*)calloc(cols,sizeof(double));
			arr[i][i] = 1;
		}
		return result;
	}


	double* operator[] (unsigned int row)
	{
		return arr[row];
	}

	
	unsigned int getRows() const
	{
		return rows;
	}
	unsigned int getCols() const
	{
		return cols;
	}

};

#endif