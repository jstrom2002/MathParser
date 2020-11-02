/**
*  \brief Parent class for parsing strings. Holds an array of declared math objects 
*	in memory as well.
*/

#pragma once
#include <string>
#include <vector>
#include "Types.h"
#include "Function.h"
#include "Matrix.h"
#include "MemoryManager.h"
#include "Polynomial.h"

namespace MathParser
{
	class Vector;
	class UserInterface;

	class MathParser
	{
	public:		
		int precision;

		MathParser();
		void clear();
		std::string parse(std::string str);

		friend class UserInterface;
	protected:
		MemoryManager memory;

	private:
		void clc();
		std::string disp(std::string str);
		std::string type(std::string str);

		// Parse input string of form "'f(x,y,z, ..)=exp(x^2)..'" into two strings: 
		// the right hand side ('exp(x^2...)...etc'), and the lhs 'x,y,z,...etc' ".
		std::vector<std::string> ParseFunction(std::string in);
		std::vector<Function> ParseNFunctions(std::string str);

		// Turns a string of format 'p1(x) = x^3 + 2x + -1' to a vector of values, 
		// then constructs a polynomial from that vector.
		Polynomial ParsePolynomial(std::string in);

		// Helper subfunctions for handling refernces to saved memory.
		std::string parseFunctionReference(std::string str, int idx);
		std::string parseMatrixReference(std::string str, int idx);
		std::string parsePolynomialReference(std::string str, int idx);
		std::string parseVectorReference(std::string str, int idx);

		// Parse matrices from the command line in the format '[1,0,0;0,1,0;0,0,1],'
		// which would be a 3x3 identity matrix.
		Matrix ParseMatrix(std::string input);
	};
}