#include "Bivector.h"

namespace MathParser
{
	Vector Bivector::operator[](int i) const
	{
		if (!i)
			return a;
		else if (i == 1)
			return b;
	}

	Vector& Bivector::operator[](int i)
	{
		if (!i)
			return a;
		else if (i == 1)
			return b;
	}
}