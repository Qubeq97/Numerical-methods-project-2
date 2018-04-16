#ifndef VECTOR
#define VECTOR
#include "Matrix.h"


class Vector :
	public Matrix
{
public:
	Vector(int);
	double& operator[] (int);
	double& operator()(int);
	double operator() (int) const;
	int getLength() const;
	Vector(const Matrix&);
};

#endif