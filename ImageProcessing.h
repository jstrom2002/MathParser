#pragma once
#include "Types.h"
#include "BaseObject.h"
#include <string>
#include <vector>

namespace MathParser
{
	class Matrix;
	class Vector;

	Matrix identityFilter(unsigned int sz);
	Matrix sharpeningFilter(unsigned int sz);
	Matrix edgeDetectorFilter(unsigned int sz);
	Matrix boxFilter(unsigned int sz);
	Matrix GaussianKernel1D(int sz, real sigma=1);
	Matrix GaussianKernel2D(int sz, real sigma=1);
	Matrix gradientFilter(std::string type, unsigned int wrt=0);
}