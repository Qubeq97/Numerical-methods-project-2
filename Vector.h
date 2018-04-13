#ifndef VECTOR
#define VECTOR
#include "Matrix.h"


class Vector :
	public Matrix
{
public:

	Vector(unsigned int elements)
	{
		this->Matrix::Matrix(elements, 1);
	}

	double& operator[] (unsigned int element)
	{
		return arr[element][0];
	}

	unsigned int getLength() const
	{
		return rows;
	}

};

#endif