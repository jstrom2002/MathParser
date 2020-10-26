/**
*  Generic base class for mathematical objects that will be saved in memory and manipulated
*  at run-time, (ie Polynomials, Functions). Shared methods are contained within 
*  this parent class.
*/

#pragma once
#include <string>

namespace MathParser
{
	class BaseObject 
	{
	public:
		
		virtual std::string to_string(int precision = 4) = 0;
		virtual real evaluate(real x, int precision = 4) { return 0; };
	};
}