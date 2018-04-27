#include "Vector.h"


Vector::Vector(int elements)
{
	this->Matrix::Matrix(elements, 1);
}



double& Vector::operator()(int element)
{
	assert(element < rows);
	return arr[element][0];
}

double Vector::operator()(int element) const
{
	assert(element < rows);
	return arr[element][0];
}

int Vector::getLength() const
{
	return rows;
}

Vector::Vector(const Matrix & m)
{
	assert(m.getCols() == 1);
	this->Matrix::Matrix(m);
}


