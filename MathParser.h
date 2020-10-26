/**
*  \brief Parent class for parsing strings. Holds an array of declared math objects in memory as well.
*/

#pragma once
#include <string>
#include <vector>
#include "Types.h"
#include "Function.h"
#include "Polynomial.h"

namespace MathParser
{
	class Function;
	class Polynomial;
	class Matrix;
	class Vector;

	class MathParser
	{
	public:		
		int precision;

		MathParser();
		std::string parse(std::string);
		std::string printSavedFunctions();
		std::string printSavedPolynomials();
		void clear();

	protected:
		real LAST;
		std::vector<char> variables;
		std::vector<Function> savedFunctions;
		std::vector<Polynomial> savedPolynomials;

	private:
		// Parse input string of form "'f(x,y,z, ..)=exp(x^2)..'" into two strings: 
		// the right hand side ('exp(x^2...)...etc'), and the lhs 'x,y,z,...etc' ".
		std::vector<std::string> ParseFunction(std::string in);
		std::vector<Function> ParseNFunctions(std::string str);

		// Turns a string of format 'p1(x) = x^3 + 2x + -1' to a vector of values, 
		// then constructs a polynomial from that vector.
		Polynomial ParsePolynomial(std::string in);

		// Helper subfunctions for handling refernces to saved memory.
		std::string parseFunctionReference(std::string str, int idx);
		std::string parsePolynomialReference(std::string str, int idx);
	};
}