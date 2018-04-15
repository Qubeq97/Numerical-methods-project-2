#ifndef VECTOR
#define VECTOR
#include "Matrix.h"


class Vector :
	public Matrix
{
public:
	Vector(unsigned int);
	double& operator[] (unsigned int);
	double& operator()(unsigned int);
	double operator() (unsigned int) const;
	unsigned int getLength() const;
	Vector(const Matrix&);
};

#endif