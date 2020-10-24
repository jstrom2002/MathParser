/**
*  List of helper functions used during string parsing for math equation parser.
*/

#pragma once
#include "Types.h"
#include <string>
#include <vector>

namespace MathParser
{
	// TESTS FOR USER-DEFINED FUNCTIONS:
	// List of all user-defined function name strings recognized by the 'evaluate()' 
	// function of the 'Function' class. This allows users to define parsing methods for
	// things like 'radix()' etc. These functions will be replaced inside the parsed string
	// during the 'Function::evalute()' call.
	extern std::vector<std::string> functionList;

	// Test a string to see if it has a user-defined function in it from the 'functionList' vector.
	int findFunction(std::string in);

	// Checks if a string has a function in it from the user-defined list of parseable functions.
	bool hasFunction(std::string in);

	// TESTS FOR OPERATORS:
	// Test a string to see if it has an operator in it (does not include parenthesis).
	bool hasOperator(std::string in);

	// Test a string to see if it has an operator in it (does not include parenthesis).
	bool hasOperatorNotComplexNumber(std::string in);

	// Gives position of nearest operator (does not include parenthesis).
	int findOperator(std::string in);

	// Returns number of operators present in a string, including parenthesis.
	int countOperators(std::string str);

	// TESTS FOR LOCATIONS OF REAL/COMPLEX NUMBERS IN A STRING:
	int findNearestReal(std::string str);
	std::string getNearestReal(std::string str);
	int countReals(std::string str);
	int findComplexNumber(std::string str);
	std::string getNearestComplexNumber(std::string str);
	std::string removeComplexNumbers(std::string str);
	std::string removeNextReal(std::string str);
	std::string removeReals(std::string str);
	int countComplexNumbers(std::string str);

	std::string convertRealToComplex(std::string funct);
	bool isOperableExpressionComplex(std::string tempstr);
	bool isOperableExpression(std::string tempstr);
	std::string processInnerTermsComplex(std::string str);
	std::string processInnerTerms(std::string str);
	std::string processTwoTermOperationComplex(std::string funct, std::string tempstr);
	std::string processTwoTermOperation(std::string funct, std::string tempstr);

	std::string getNextInnerMostExpression(std::string str);
	std::vector<std::string> loadExpressionsIntoArray(std::string str);

	// GENERAL FUNCTIONS TO CONVERT STRINGS OF ARITHMETIC STATEMENTS TO NUMBERS:
	complex calculateArithmeticComplex(std::string funct);
	real calculateArithmetic(std::string funct);

	// FUNCTIONS FOR PARSING NUMBERS:
	// Wrapper functions around C++'s string to real parsing functions, ie 'std::stod(),'
	// to format input and prevent exceptions.
	real ParseInputNumber(std::string in);
	complex ParseComplexNumber(std::string in);

	// Parses a comma-delimited array of real/complex numbers.
	std::vector<real> ParseNInputNumbers(std::string str);	
	std::vector<complex> ParseNComplexNumbers(std::string str);
}