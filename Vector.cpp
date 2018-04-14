#include "Vector.h"


Vector::Vector(unsigned int elements)
{
	this->Matrix::Matrix(elements, 1);
}

double& Vector::operator[] (unsigned int element)
{
	return arr[element][0];
}

unsigned int Vector::getLength() const
{
	return rows;
}

Vector::Vector(const Matrix & m)
{
	assert(m.getCols() == 1);
	this->Matrix::Matrix(m);
}


