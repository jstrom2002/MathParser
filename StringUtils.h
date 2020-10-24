/**
*   Utility functions for use with math parser.
*/

#pragma once
#include <string>
#include <vector>
#include "Types.h"

namespace MathParser
{
	bool isOperator(char c);
	bool isOperatorNotParen(char c);
	char randomChar();

	std::string eraseParenPair(std::string str, size_t pos);
	size_t findMatchingParen(std::string str, int pos);
	std::string getMantissa(std::string n);
	std::string getInteger(std::string n);
	std::string getNextInnerMostExpression(std::string str);
	bool isInteger(std::string n);
	std::vector<std::string> loadExpressionsIntoArray(std::string str);
	std::string matchParens(std::string str);
	std::string removeUnmatchedParens(std::string funct);
	std::string replaceChar(std::string str, char replace, char replacer);
	std::string replaceString(std::string str, std::string replace, std::string replacer);
	std::string removeSpaces(std::string str);
	std::string reverseString(std::string str);
	std::string to_stringPrecision(real val, unsigned int precisn);
	std::string toHex(int val);
	std::vector<std::string> tokenize(std::string toTokenize, std::string token);
	std::string trim(std::string str);
}