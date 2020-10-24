/**
*  Helper functions for dealing with complex numbers.
*/

#pragma once
#include "Types.h"

namespace MathParser
{
	std::string toString(complex z, int precision = 4);
	std::string toStringBothParts(complex z, int precision = 4);
}