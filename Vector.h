#ifndef VECTOR
#define VECTOR
#include "Matrix.h"


class Vector :
	public Matrix
{
public:
	Vector(unsigned int elements);
	double& operator[] (unsigned int element);
	unsigned int getLength() const;
	Vector(const Matrix&);
};

#endif