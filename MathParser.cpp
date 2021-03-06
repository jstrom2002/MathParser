#include "MathParser.h"
#include "StringUtils.h"
#include "Polynomial.h"
#include "ComplexNumber.h"
#include "Function.h"
#include "Parsing.h"
#include "StringUtils.h"
#include "Calculus.h"
#include "FileIO.h"
#include <string>
#include <algorithm>
#include <regex>
#include <iostream>
#include <cstdlib>
#include <fstream>

namespace MathParser
{
	MathParser::MathParser() : precision(4)
	{
	}

	void MathParser::clc()
	{
#if defined(WINDOWS) || defined(_WIN32) || defined(_WIN64)
		std::system("cls");
#else
		// Assume POSIX
		std::system("clear");
#endif
	}

	std::string MathParser::disp(std::string str)
	{
		// Remove ')' around outside of filename.
		while (str[str.length() - 1] == ')')
			str = str.substr(0, str.length() - 1);

		// Remove all chars before "(" from from beginning of filename.
		str = str.substr(str.find("(") + 1);

		// Remove all spaces from inner string.
		str = removeSpaces(str);

		// First char will always represent a type of object, ie 'M0' is always
		// a matrix.
		char type = str[0];
		int idx = (int)ParseInputNumber(str.substr(1));
		if (idx < 0)
			return "";

		switch (type)
		{
		case 'f'://Function
			if (idx < memory.functions.size())
				return memory.functions[idx].to_string(this->precision);
		case 'M'://Matrix
			if (idx < memory.matrices.size())
				return "\n" + memory.matrices[idx].to_string(this->precision);				
		case 'p'://Polynomial
			if (idx < memory.polynomials.size())
				return memory.polynomials[idx].to_string(this->precision);
		case 'v'://Vector
			if (idx < memory.vectors.size())
				return memory.vectors[idx].to_string(this->precision);
		}

		return "";
	}

	std::string MathParser::parseFunctionReference(std::string str, int idx)
	{
		if (idx >= memory.functions.size())
			return "";

		if (str.find("evaluate(") != std::string::npos)
		{
			// Get terms in parentheses and evaluate.
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			memory.LAST = memory.functions[idx].evaluate(vals, this->precision);
			return to_string_precision(memory.LAST, this->precision);
		}
		else if (str.find("derivative(") != std::string::npos)
		{
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			memory.LAST = nthDerivative1D(vals[0],vals[1],memory.functions[idx], vals[2]);
			return to_string_precision(memory.LAST, this->precision);
		}
		else if (str.find("integral(") != std::string::npos)
		{
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			//LAST = integrateSimpsonsRule(vals[0], vals[1], savedFunctions[idx], vals[2]);
			//LAST = integrateGaussLegendreQuadrature(vals[0], vals[1], savedFunctions[idx]);
			memory.LAST = integrateGaussKronodQuadrature(vals[0], vals[1], memory.functions[idx]);
			return to_string_precision(memory.LAST, this->precision);
		}

		return "";
	}

	std::string MathParser::parseMatrixReference(std::string str, int idx)
	{
		if (idx >= memory.matrices.size())
			return "";

		if (str.find("reshape(") != std::string::npos)
		{
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			memory.matrices[idx].reshape(vals[0],vals[1]);
		}
		else if (str.find("submatrix(") != std::string::npos)
		{
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			memory.matrices.push_back(memory.matrices[idx].submatrix(vals[0], vals[1],
				vals[2],vals[3]));
		}
		else if (str.find("inverse()") != std::string::npos)
		{
			memory.matrices.push_back(memory.matrices.back().inverse());
			return memory.matrices.back().to_string(precision);
		}
		else if (str.find("det()") != std::string::npos)
		{
			return to_string_precision(memory.matrices.back().det(), precision);
		}

		return "";
	}

	std::string MathParser::parsePolynomialReference(std::string str, int idx)
	{
		if (idx >= memory.polynomials.size())
			return "";

		if (str.find("evaluate(") != std::string::npos)
		{
			// Get terms in parentheses and evaluate
			std::vector<real> vals = ParseNInputNumbers(
				str.substr(str.find("("), str.find(")")));
			memory.LAST = memory.polynomials[idx].evaluate(vals[0], this->precision);
			return to_string_precision(memory.LAST, this->precision);
		}
		else if (str.find("realRoots()") != std::string::npos)
		{			
			std::vector<real> vals = memory.polynomials[idx].realRoots();
			std::string outputStr = "real roots: ";
			for (int i = 0; i < vals.size(); ++i)
			{
				outputStr.append(to_string_precision(vals[i], precision));
				if (i < vals.size() - 1)
					outputStr.append(",");
			}
			return outputStr;
		}
		else if (str.find("roots()") != std::string::npos)
		{
			std::vector<ComplexNumber> vals = memory.polynomials[idx].roots();
			std::string outputStr = "roots: ";
			for (int i = 0; i < vals.size(); ++i)
			{
				outputStr.append(to_string_precision(vals[i], precision));
				if (i < vals.size() - 1)
					outputStr.append(",");
			}
			return outputStr;
		}

		return "";
	}

	std::string MathParser::parseVectorReference(std::string str, int idx)
	{
		if (idx >= memory.vectors.size())
			return "";

	}

	void MathParser::clear()
	{
		memory.clear();
	}	

	std::string MathParser::parse(std::string str)
	{
		if (!str.length() || str == " " || str == "\n") 
			return ""; 

		// Replace all occurances of word 'LAST' with value for 'last.'
		if (str.find("LAST") != std::string::npos)
		{
			std::string laststr = to_string_precision(memory.LAST, this->precision);
			str = replaceString(str, "LAST", laststr);
		}

		// Check for reference to saved objects, referenced by a letter then a number.
		// (ie 'M0.eigenvalues()').
		if (std::regex_match(str, std::regex("[a-z][0-9]\..*")) ||
			std::regex_match(str, std::regex("[a-z][0-9][0-9]\..*")) ||
			std::regex_match(str, std::regex("[A-Z][0-9]\..*")) ||
			std::regex_match(str, std::regex("[A-Z][0-9][0-9]\..*")))
		{
			// Get index of referenced object.
			std::string prestr = str.substr(0, str.find(")"));
			int idx = std::stoi(str.substr(1, str.find("(")));
			if (idx < 0)
				return "";
			std::string afterDot = str.substr(str.find(".") + 1);

			// Find reference to saved objects.
			if (prestr[0] == 'f')
				return parseFunctionReference(afterDot, idx);
			else if (prestr[0] == 'p')
				return parsePolynomialReference(afterDot, idx);
			else if (prestr[0] == 'M')
				return parseMatrixReference(afterDot, idx);
			else if (prestr[0] == 'v')
				return parseVectorReference(afterDot, idx);
		}

		// Check for declared functions, polynomials, etc.
		else if (str.find("=") != std::string::npos)
		{
			// Save function or polynomial.
			if (str.find("p(") != std::string::npos)			
				memory.polynomials.push_back(ParsePolynomial(str));			
			else			
				memory.functions.push_back(ParseFunction(str));	

			return "";
		}

		// Check for declared matrices.
		else if (str.find("[") != std::string::npos && str.find("]") != std::string::npos)
		{
			memory.matrices.push_back(ParseMatrix(str));
			return "";
		}

		// If no functions, then parse arithmetic string.
		if (hasOperator(str)) 
		{
			// Remove all spaces in input to format input properly.
			str = removeSpaces(str);

			Function fp("", str);
			memory.LAST = fp.evaluate(0, this->precision);
			return to_string_precision(memory.LAST, this->precision);
		}
		else
			return "";
	}

	std::vector<std::string> MathParser::ParseFunction(std::string in) 
	{
		std::string str1;
		std::string str2;
		std::vector<std::string> answer;

		int leftbracket = 0;
		int rightbracket = 0;
		int EQsign = 0;
		bool needLB = true;
		bool needRB = true;
		bool needEQ = true;
		for (int i = 0; i < in.length(); ++i) {
			if (in[i] == ')' && needRB) 
			{
				rightbracket = i;
				needRB = false;
			}
			if (in[i] == '(' && needLB) 
			{
				leftbracket = i;
				needLB = false;
			}
			if (in[i] == '=' && needEQ) 
			{
				EQsign = i;
				needEQ = false;
			}
		}
		if (needLB || needRB || needEQ) 
			return answer;

		str1 = in.substr(EQsign + 1);//actual function
		str2 = in.substr(leftbracket + 1, rightbracket - 2);//list of variables
		
		// Trim empty spaces from strings.
		str1 = trim(str1);
		str2 = trim(str2);

		//Function is declared like so:  f[0] = function, f[1] = variables.
		answer.push_back(str2);
		answer.push_back(str1);
		return answer;
	}

	Polynomial MathParser::ParsePolynomial(std::string input) 
	{
		input = removeSpaces(input);

		// Cut off anything proceeding the 'x^3 + ...' part like a 'p1(x) = ' statement.
		if (input.find("=") != std::string::npos) 
			input = input.substr(input.find("=") + 1);		
		
		// Remove blank space at the beginning of the string
		while (input[0] == ' ')
			input = input.substr(1);

		// Insert a '+' in front of any '-' sign in case there isnt' one already there
		for (int i = 0; i < input.length(); ++i) {
			if (input[i] == '-') {
				
				if (std::isalpha(input[i + 1])) 
					input.insert(input.begin() + i + 1, '1');
				
				//go backwards from i and check if you hit a number.  if so, add a '+'
				for (int k = i; k >= 0; --k) 
				{
					if (std::isdigit(input[k]))
					{
						input.insert(input.begin() + i, '+');
						break;
					}
					if (input[k] == '+') 
						break;
				}
			}
		}

		// Find terms where exponent >= 2.
		std::vector<real> coefs;
		std::vector<int> expos;
		while (input.find("+"))
		{
			// Get everything up to next '+' symbol, which even negative numbers
			// will have now.
			size_t nextIdx = std::string::npos;
			std::string tmp = input.substr(0, nextIdx=input.find("+"));
			
			if (nextIdx == std::string::npos || !tmp.length())
				break;

			// Catch case where polynomial term is scalar.
			if (tmp.find("x") == std::string::npos)
			{
				coefs.push_back(ParseInputNumber(tmp));
				expos.push_back(0);
			}

			// Catch case where there is no exponent (ie '3x').
			else if (tmp.find("^") == std::string::npos)
			{
				std::string coef = tmp.substr(0, tmp.find("x"));
				coefs.push_back(coef.length() ? ParseInputNumber(coef) : 1.0);
				expos.push_back(1);
			}

			// Else, term should be formatted as 'ax^b'.
			else
			{
				std::string coef = tmp.substr(0, tmp.find("x"));
				coefs.push_back(coef.length() ? ParseInputNumber(coef) : 1.0);
				std::string exp = tmp.substr(tmp.find("^") + 1);
				expos.push_back(exp.length() ? ParseInputNumber(exp) : 1.0);
			}

			// Move down string to next term.
			if(nextIdx != std::string::npos)
				input = input.substr(nextIdx + 1);
		}

		// Catch case where string still has a final, scalar polynomial term.
		if (input.length())
		{
			coefs.push_back(ParseInputNumber(input));
			expos.push_back(0);
		}

		// Get largest value in exponent array to determine max # coefficients.
		int maxExponent = *std::max_element(expos.begin(), expos.end());
		if (maxExponent > 0)
			maxExponent++;

		// Create an array of coefficients indexed by exponent.
		std::vector<real> v;
		v.resize(maxExponent,0);	
		for (int i = 0; i < expos.size(); ++i)
		{
			v[expos[i]] = coefs[i];
		}
				
		// Finally, reverse array (since 'x^0' will be the first element in the
		// array, yet the first value in the string will be the highest exponent).
		std::reverse(v.begin(), v.end());

 		Polynomial p(v);
		return p;
	}

	std::vector<Function> MathParser::ParseNFunctions(std::string str) {
		str = removeSpaces(str);
		std::vector<Function> answer;
		std::vector<int> position;

		//get the number of values in the string by finding all ',' chars
		for (int i = 0; i < str.length(); ++i) {
			if (str[i] == ',') {
				position.push_back(i);
			}
			if (str[i] == ')') {//break if you've hit the end
				i = str.length();
			}
		}
		if (!position.size())
		{ //if there is only one single value in the string
			answer.push_back(Function(ParseFunction(str))); 
			return answer; 
		}
		if (position.size() > 0) {//if there are multiple values in the string
			std::string tempStr = str.substr(0, position[0]);   //get first substring
			answer.push_back(Function(ParseFunction(tempStr)));
			for (int i = 1; i < position.size(); ++i) {
				tempStr = str.substr(position[i - 1] + 1, position[i] - position[i - 1] - 1);
				answer.push_back(Function(ParseFunction(tempStr)));
			}
			tempStr = str.substr(position[position.size() - 1] + 1);
			answer.push_back(Function(ParseFunction(tempStr)));
		}
		return answer;
	}

	Matrix MathParser::ParseMatrix(std::string input)
	{
		std::vector<real> vals;
		int rows = 0;
		int cols = 0;

		input = removeSpaces(input);

		// Remove beginning and ending brackets (']' chars).
		while (input[0] == '[')
			input = input.substr(1);
		while (input[input.length()-1] == ']')
			input = input.substr(0, input.length()-1);

		// Tokenize by end of row marker ';'.
		std::vector<std::string> strs = tokenize(input, ";");
		rows = strs.size();
		for (int i = 0; i < strs.size(); ++i)
		{
			// Parse the next row's string, then add those values to the output array.
			std::vector<real> temp = ParseNInputNumbers(strs[i]);
			if (temp.size() > cols)
				cols = temp.size();
			vals.insert(std::end(vals), std::begin(temp), std::end(temp));
		}

		return Matrix(rows, cols, vals);
	}

	std::string MathParser::type(std::string str)
	{
		str = str.substr(str.find(" ")+1);

		std::ifstream filep(str, std::ios::out);
		if (!filep.is_open())
			return "";

		std::string outputStr = "";
		std::string buff = "";
		while (getline(filep, buff))
			outputStr.append(buff);

		return outputStr;
	}
}