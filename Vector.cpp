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