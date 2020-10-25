/**
*   Class for parsing/handling user-input functions. By default, even simple mathematical expressions
*	are parsed by members of this class.
*/

#pragma once
#include "Types.h"
#include <vector>
#include "BaseObject.h"

namespace MathParser
{
	class Function : public BaseObject
	{
	public:
		Function(std::string str);
		Function(std::string vars, std::string funct);
		Function(std::vector<char> vars, std::string funct);
		Function(std::vector<std::string> funct);

		real evaluate(std::vector<real> vals, int precision = 4);
		virtual real evaluate(real x, int precision = 4);
		virtual std::string to_string(int precision = 4);

	private:

		// Prevent default construction method.
		Function() = default;

		// Strictly contains variable chars from lhs of function, i.e. it will contain 
		// 'x,y' for function 'f(x,y) = x^2 + y^2'.
		std::vector <char> variables;

		// Contains just the rhs of the the function (using the example above, it would be 'x^2 + y^2').
		std::string function;

		// Convert string to real if the string is a function of the kind that is parsable.
		real evalFunction(std::string str);

		// Helper function for parsing input variables of a function to an array of chars.
		std::vector<char> varsToChars(std::string vars);
	};
}