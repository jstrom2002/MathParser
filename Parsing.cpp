#include "Parsing.h"
#include "StringUtils.h"
#include "MathLib.h"
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

namespace MathParser
{
	// Declare 'extern' object here for program-wide use.
	std::vector<std::string> functionList = 
	{
		"abs", "acos", "acosh", "arg", "asin", "asinh", "atan", "atan2", "atanh", 
		"BesselFunctionFirstKind", "BesselFunctionSecondKind", "beta", "cbrt", "Ci", "conj", 
		"ceil", "combination", "cos", "cosh", "cot", "csc", "digamma", "doubleFactorial", 
		"Ei", "erf", "erfc", "exp", "exp2", "factorial", "fallingFactorial", "floor", 
		 "gamma", "hypergeometricFunction", "hypot", "if", "imag", "int", "inverseErfc", 
		 "LambertW", "Li", "log", "log2", "log10", "logGamma", "logit", "logisticFunction",
		"max", "min", "permutation", "prime", "sumPrimes", "productPrimes", "primeCountingFunction",
		"polar", "pow", "polygamma", "real", "regularModifiedBesselFunctionFirstKind",
		"rect", "risingFactorial", "sec", "Si", "sin", "sinc", "sinh", 
		"sphericalBesselFunctionFirstKind", "squarewave", "sgn", "sqrt", "step", 
		"StirlingNumberFirstKind", "StirlingNumberSecondKind", "tan", "tanh", "tri", 
		"trunc", "WrightOmegaFunction", "zeta"
	};

	int countOperators(std::string str)
	{
		// Format string by adding '+' signs to the left of negative signs. In all cases,
		// the negative sign should be treated as a modifier to a number, not an operator.
		// While '+' is only an operator. So doing this preserves the intent of the original
		// mathematical statement.
		int nearestNegative = std::string::npos;
		while ((nearestNegative = str.find("-")) > 0 && str[nearestNegative-1] != '+' )
		{
			str.insert(str.begin()+nearestNegative, '+');
		}

		if (findComplexNumber(str) > 0) 
			str = removeComplexNumbers(str);
		
		str = removeReals(str);
		
		if (str == "" || str.size() == 0) //handle null string
			return 0;

		int count = 0;
		for (int i = 0; i < str.length(); ++i) 
		{
			if (isOperatorNotParen(str[i])) 
				++count;
		}

		return count;
	}

	bool hasOperatorNotComplexNumber(std::string in)
	{
		if (in.size() == 0) 
			return false;

		in = removeComplexNumbers(in);

		if (//find defined operators
			in.find("*") != std::string::npos ||
			in.find("/") != std::string::npos ||
			in.find("^") != std::string::npos ||
			in.find("%") != std::string::npos
			) {
			return true;
		}

		//defined functions
		for (int i = 0; i < functionList.size(); ++i) {
			if (in.find(functionList[i]) != std::string::npos) {
				return true;
			}
		}

		//finally, check to see if there is an operator not attached to a complexNumber
		if (in.find("+") != std::string::npos) {//check '+' operator, rejecting any complex number of form '1+2i'
			int i = in.find('+');
			for (i + 1; i < in.length(); ++i) {
				if (std::isdigit(in[i]) == false && in[i] != ',') {
					if (in[i] == 'i') {
						in = in.substr(i);
						i = in.find('+');
					}
					else {
						return true;
					}
				}
			}
		}

		if (in.find("-") != std::string::npos) 
		{
			bool halfway = false;
			int i = in.find('-');
			for (i + 1; i < in.length(); ++i) 
			{
				if (halfway && i == in.length() - 1)
					return true;
				if (!std::isdigit(in[i]) && in[i] != '.' && in[i] != ',') {
					if (!halfway && (in[i] == '-' || in[i] == '+')) 
						halfway = true;
					if (in[i] == 'i') 
					{
						in = in.substr(i + 1);
						i = in.find('+');
						int i2 = in.find('-');
						if (i2 >= 0 && i2 < i) { i = i2; }
						halfway = false;
					}
				}
			}
		}

		return false;	//if nothing is thrown, return false by default
	}

	int findOperator(std::string in)
	{
		int a = std::string::npos;
		if (in.find("+") >= 0) { a = ((int)in.find("+") < a) ? in.find("+") : a; }
		if (in.find("-") >= 0) { a = ((int)in.find("-") < a) ? in.find("-") : a; }
		if (in.find("*") >= 0) { a = ((int)in.find("*") < a) ? in.find("*") : a; }
		if (in.find("/") >= 0) { a = ((int)in.find("/") < a) ? in.find("/") : a; }
		if (in.find("^") >= 0) { a = ((int)in.find("^") < a) ? in.find("^") : a; }
		if (in.find("%") >= 0) { a = ((int)in.find("%") < a) ? in.find("%") : a; }
		return a;
	}

	bool hasOperator(std::string in) 
	{
		if (
			//defined operators
			in.find("+") != std::string::npos ||
			in.find("-") != std::string::npos ||
			in.find("*") != std::string::npos ||
			in.find("/") != std::string::npos ||
			in.find("^") != std::string::npos ||
			in.find("%") != std::string::npos
			) {
			return true;
		}

		// Check for user-defined functions
		for (int i = 0; i < functionList.size(); ++i) {
			if (in.find(functionList[i]) != std::string::npos) {
				return true;
			}
		}

		// If nothing is thrown, return false by default.
		return false;	
	}


	int findFunction(std::string in) 
	{
		for (int i = 0; i < functionList.size(); ++i) {
			if (in.find(functionList[i]) != std::string::npos) {
				return i;//return index of function in array of pre-defined functions
			}
		}
		return -1;	//if nothing is thrown, return -1 by default
	}

	bool hasFunction(std::string in) 
	{
		for (int i = 0; i < functionList.size(); ++i) {
			if (in.find(functionList[i]) != std::string::npos) {
				return true;
			}
		}
		return false;	//if nothing is thrown, return false by default
	}

	int findNearestReal(std::string str) 
	{
		// Format string to function properly.
		str = removeSpaces(str);

		int start = 0;
		for (int i = 0; i < str.length(); ++i) 
		{
			if (std::isdigit(str[i])) 
			{
				bool isNegativeNumber = false;
				if (i > 0 && str[i - 1] == '-') 
				{
					int j = i - 1;
					while (j >= 0 && str[j] != '(' && !std::isdigit(str[j])) 
					{
						if ((j > 0 && str[j - 1] == '(') || j == 0) 
						{
							isNegativeNumber = true; 
							j = -1; 
						}
						--j;
					}
				}
				start = i;
				if (isNegativeNumber) 
					start = i - 1;
				return start;
			}
		}
		return std::string::npos;
	}

	// Wrapper around std::stod to allow user to change parsing method.
	real ParseInputNumber(std::string in) 
	{ 
		// Format string.
		for (int i = 0; i < in.length(); ++i)
		{
			if (!std::isdigit(in[i]) && in[i] != '-' && in[i] != '.')
			{
				in.erase(in.begin() + i);
				i--;
			}
		}

		if (!in.length())
			return 0;

		return std::stod(in); 
	}

	std::vector<complex> ParseNComplexNumbers(std::string str) {
		str = removeSpaces(str);
		std::vector<complex> answer;
		std::vector<std::string> strs = tokenize(str, ",");
		for (int i = 0; i < strs.size(); ++i) {
			if (findComplexNumber(strs[i]) != std::string::npos) {
				answer.push_back(ParseComplexNumber(strs[i]));
			}
			else if (findNearestReal(strs[i]) != std::string::npos) {
				answer.push_back(complex(ParseInputNumber(strs[i])));
			}
		}
		return answer;
	}

	int findComplexNumber(std::string str) {//return index if ComplexNumber is found, -1 if not
		
		int idx = findNearestReal(str);

		// A complex number needs more than one digit in it, so if it is a single value return -1.
		// This method does not allow 'i' itself to be a complex number.
		if (idx >= str.length() - 1)
			return -1;

		while (idx >= 0)
		{
			// while next value is a digit (or special char), keep searching.
			int i = idx + 1;
			while (i < str.size() && (std::isdigit(str[i]) || str[i] == '.' || 
				str[i] == '-' || str[i] == '+' || str[i] == ' '))
			{
				i++;
			}

			// If final 'i' value is found, return.
			if (str[i] == 'i')
			{
				return idx;
			}

			// Else, keep going.
			else
			{
				str = str.substr(i);
				idx = findNearestReal(str);
			}
		}

		return -1;
	}

	std::string getNearestComplexNumber(std::string str) {
		int i = findComplexNumber(str);
		int i2 = str.substr(i).find("i") + 1;
		if (str.substr(i).find("i") < 0) { return NULL; }
		else {
			std::string cn = str.substr(i, i2);
			while (cn[0] != '-' && std::isdigit(cn[0]) == false) { cn = cn.substr(1); }
			while (cn[cn.length() - 1] != 'i') { cn.pop_back(); }
			return cn;
		}
		return NULL;
	}

	std::string getNearestReal(std::string str) 
	{
		if (!str.length())
			return "";

		// Format string to work properly.
		str = removeSpaces(str);

		// Catch case where number is first and only digit, or last digit.
		int i = findNearestReal(str);
		if (str.length() == 1 && i == 0)
		{
			std::string temp = " ";
			temp[0] = str[0];
			return str;
		}
		else if (i == str.length() - 1)
		{
			if (str[i - 1] == '-')
				return str.substr(i - 1);
			else
				return str.substr(i);
		}

		// Search string for last digit.
		std::string returnStr = " ";
		returnStr[0] = str[i];
		i++;
		while (i < str.size() && (std::isdigit(str[i]) || str[i] == '.'))
		{
			returnStr.push_back(str[i]);
			i++;
		}

		// Remove values that might be incorrect from the string.
		while (returnStr[0] != '-' && !std::isdigit(returnStr[0]))
			returnStr = returnStr.substr(1);
		while (returnStr[returnStr.length() - 1] && !std::isdigit(returnStr[returnStr.length() - 1]))
			returnStr = returnStr.substr(0, returnStr.length() - 1);

		return returnStr;
	}

	std::string removeNextReal(std::string str)
	{
		int idx = -1;
		if ((idx = findNearestReal(str)) < 0)
			return str;

		std::string cpnum = getNearestReal(str);
		str.erase(idx, cpnum.length());		
		return str;
	}

	std::string removeReals(std::string str) {

		int idx = -1;
		if ( (idx = findNearestReal(str)) < 0) 
			return str;
		
		while (idx >= 0) 
		{
			std::string cpnum = getNearestReal(str);
			str.erase(idx, cpnum.length());
			idx = findNearestReal(str);
		}

		return str;
	}

	std::string removeComplexNumbers(std::string str) 
	{
		if (findComplexNumber(str) < 0) 
			return str;
		
		while (findComplexNumber(str) >= 0) 
		{
			std::string cpnum = getNearestComplexNumber(str);
			str.erase(findComplexNumber(str), cpnum.length());
		}

		return str;
	}

	int countReals(std::string str) {//returns number of ComplexNumbers present in a string
		if (str == "" || str.size() == 0) //handle null string
			return 0; 

		int count = 0;
		while (findNearestReal(str) >= 0) 
		{
			std::string cpnum = getNearestReal(str);
			size_t end = cpnum.size();
			size_t start = findNearestReal(str);
			str.erase(findNearestReal(str), cpnum.size());
			++count;
		}
		return count;
	}

	int countComplexNumbers(std::string str) {//returns number of ComplexNumbers present in a string
		if (str == "" || str.size() == 0) //handle null string
			return 0; 

		// Format string to function properly.
		str = removeSpaces(str);

		int count = 0;
		while (findComplexNumber(str) >= 0) 
		{
			std::string cpnum = getNearestComplexNumber(str);
			str.erase(findComplexNumber(str), cpnum.size());
			++count;
		}

		return count;
	}

	std::string convertRealToComplexNumber(std::string str, int start, int end) 
	{
		std::string replacer = "(";
		complex z = complex(ParseInputNumber(str.substr(start, end - start)));
		replacer.append(to_string_precision(z, 4));
		replacer.append(")");
		str.erase(start, end - start);
		str.insert(start, replacer);
		return str;
	}

	std::string convertRealToComplex(std::string funct) 
	{
		//first, convert any real values in the string to complex
		bool foundNum = false;
		int start = 0;
		int end = 0;
		for (int i = 0; i < funct.length(); ++i) {
			if (funct[i] == 'i') {	//if it is a complex number, skip
				start = end = i;
				foundNum = false;
			}
			if (i == funct.length() - 1 && foundNum && std::isdigit(funct[i])) {
				end = i;
				funct = convertRealToComplexNumber(funct, start, end);
				start = end = i;
				foundNum = false;
			}
			if ( !std::isdigit(funct[i]) && start != end && funct[i] != '.' 
				&& findComplexNumber(funct.substr(i - 1)) != 0) {//find end of real
				if (foundNum) {
					end = i;
					funct = convertRealToComplexNumber(funct, start, end);
				}
				start = end = i;
				foundNum = false;
			}
			if (std::isdigit(funct[i])) {//if the char is a number, include it in the range
				bool isNegativeNumber = false;
				if (i > 0 && funct[i - 1] == '-') {
					int j = i - 1;
					while (j > 0 && funct[j] != '(' && !std::isdigit(funct[j])) {
						if (funct[j - 1] == '(') { isNegativeNumber = true; j = -1; }
						--j;
					}
				}

				if (start == end) {
					start = i;
					if (isNegativeNumber) 					
						start = i - 1;					
					foundNum = true;
				}
				end = i + 1;
				if (i == funct.length() - 1 && foundNum) {
					funct = convertRealToComplexNumber(funct, start, end);
					foundNum = false;
					i += end - start;
				}
			}
		}
		return funct;
	}

	std::string processTwoTermOperation(std::string funct, std::string tempstr) {
		real x = calculateArithmetic(tempstr);
		std::string replc = ("(");
		replc.append(std::to_string(x));
		replc.append(")");
		funct = replaceString(funct, tempstr, replc);
		return funct;
	}

	std::string processTwoTermOperationComplex(std::string funct, std::string tempstr) {
		complex x = calculateArithmeticComplex(tempstr);
		std::string replc = ("(");
		replc.append(to_string_precision(x, 4));
		replc.append(")");
		funct = replaceString(funct, tempstr, replc);
		return funct;
	}

	bool isOperableExpression(std::string tempstr) 
	{
		if (countReals(tempstr) == 2 && countOperators(tempstr) == 1)
			return true;	

		return false;
	}

	bool isOperableExpressionComplex(std::string tempstr) {
		if (countComplexNumbers(tempstr) == 2 && countOperators(tempstr) == 1) 		
			return true; 		

		return false;
	}

	std::string processInnerTerms(std::string str) {
		if (str.find('(') == std::string::npos) { return str; }
		bool hasOperableTerms = true;
		while (hasOperableTerms) 
		{
			std::vector<std::string> exps = loadExpressionsIntoArray(str);
			for (int i = 0; i < exps.size(); ++i) {
				if (countReals(exps[i]) >= 2 && countOperators(exps[i]) >= 1) 
				{
					str = processTwoTermOperation(str, exps[i]);
					i = exps.size() + 1;
				}
				if (i == exps.size() - 1)
					hasOperableTerms = false;
			}
		}
		return str;
	}

	std::string processInnerTermsComplex(std::string str) {
		bool hasOperableTerms = true;
		while (hasOperableTerms) 
		{
			std::vector<std::string> exps = loadExpressionsIntoArray(str);
			for (int i = 0; i < exps.size(); ++i) {
				if (countComplexNumbers(exps[i]) >= 2 && countOperators(exps[i]) >= 1) {
					str = processTwoTermOperationComplex(str, exps[i]);
					i = exps.size() + 1;
				}
				if (i == exps.size() - 1) 
					hasOperableTerms = false;
			}
		}
		return str;
	}


	std::string getNextInnerMostExpression(std::string str) {
		/* Helper function for parser. This will get all parens and load their locations
		into vectors for comparison, and also assign values to each closed paren of
		how many outer brackets are enclosed within it (ie the first term of
		'(((x-2) +...))' has a value of 3 outer brackets
		ex: "(((1-2)*(3-1)-(10))-((sin(1)*2)))"
		openParen = {0,1,2,8,14,20,21,25}
		closeParen = {6,12,17,18,27,30,31,32}
		values = {3,3,3,2,4,3,2,1}
		largest value = 4, start = will return "(1)"	*/
		std::vector<size_t> openParen, closeParen;
		std::vector<int> values;
		int val = 0;
		for (int i = 0; i < str.length(); ++i) {
			if (str[i] == '(') { ++val; openParen.push_back(i); }
			if (str[i] == ')') { values.push_back(val); --val; closeParen.push_back(i); }
		}

		int maxParenValIndex = std::max_element(values.begin(), values.end()) - values.begin();
		if (maxParenValIndex > 0) {
			size_t end = closeParen[maxParenValIndex];
			size_t start = findMatchingParen(str, end);
			return str.substr(start, end - start + 1);
		}
		else { return NULL; }
	}

	std::vector<std::string> loadExpressionsIntoArray(std::string str) {
		/* Helper function for parser. This will get all parens and load
		their locations into vectors for comparison, and also
		assign values to each closed paren of how many outer brackets are enclosed within it
		(ie the first term of'(((x-2) +...))' has a value of 3 outer brackets
		ex: "(((1-2)*(3-1)-(10))-((sin(1)*2)))"
		openParen = {0,1,2,8,14,20,21,25}
		closeParen = {6,12,17,18,27,30,31,32}
		values = {3,3,3,2,4,3,2,1}
		largest value = 4, will return 8 parenthetical expressions (since
		there are 8 paren pairs)	*/
		std::vector<int> openParen, closeParen;
		std::vector<int> values;
		std::vector<std::string> expression;
		int val = 0;
		for (int i = 0; i < str.length(); ++i) {
			if (str[i] == '(') { ++val; openParen.push_back(i); }
			if (str[i] == ')') { values.push_back(val); --val; closeParen.push_back(i); }
		}

		for (int i = 0; i < closeParen.size(); ++i) {
			int maxParenValIndex = std::max_element(values.begin(), values.end()) - values.begin();
			if (maxParenValIndex >= 0 && maxParenValIndex < closeParen.size()) {
				size_t end = closeParen[maxParenValIndex];
				size_t start = findMatchingParen(str, end);
				expression.push_back(str.substr(start, end - start + 1));
				closeParen[maxParenValIndex] = -9999;
				values[maxParenValIndex] = -9999;
			}
		}
		return expression;
	}

	complex calculateArithmeticComplex(std::string funct) {
		int start = 0;
		int end = 0;
		std::vector<complex> dbls;
		std::vector<char> ops;

		//load all complex numbers into the array
		while (findComplexNumber(funct) >= 0) {
			std::string cpnum = getNearestComplexNumber(funct);
			dbls.push_back(complex(ParseComplexNumber(cpnum)));
			funct.erase(findComplexNumber(funct), cpnum.size());
		}
		if (dbls.size() == 1) { return dbls[0]; }

		//collect leftover operators for arithmetic evaluation
		while (hasOperator(funct)) {
			for (int i = 0; i < funct.length(); ++i) {
				if (isOperatorNotParen(funct[i])) 
				{
					ops.push_back(funct[i]);
					funct.erase(i, 1);
				}
			}
		}
		funct.clear();

		//if operator vector is empty, there are no operators, only a value to return
		if (ops.size() == 0) { return ParseComplexNumber(funct); }

		complex answer(0);
		int pass = 0;

		//evaluate in PEMDMAS order, 6 passes through the vector, eliminating evaluated terms and operators as you go
		while (dbls.size() > 1) {

			for (int i = 0; i < ops.size(); ++i) {
				//note: format = 'a+b'
				complex a = dbls[i];
				complex b = dbls[i + 1];
				if (i >= 0 && ops[i] == '+' && pass == 4) {
					answer = a + b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '-' && pass == 5) {
					answer = a - b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '*' && pass == 1) {
					answer = a * b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '/' && pass == 2) {
					if (b == complex(0)) {
						return (complex());
					}
					else {
						answer = a / b;
					}
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '^' && pass == 0) {
					answer = pow(a, b);
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				//if (i >= 0 && ops[i] == '%' && pass == 3) {
				//	answer = (int)a.modulus() % (int)b.modulus();
				//	dbls[i + 1] = answer;
				//	dbls.erase(dbls.begin() + i);
				//	ops.erase(ops.begin() + i);
				//	--i;
				//}
			}
			++pass;
			pass %= 6;//keep going around until only one dbl remains
		}
		return answer; //note: if unsuccessful, will return 0
	}

	real calculateArithmetic(std::string funct) {

		//// Format input string to prevent issues during calculation.
		funct = removeSpaces(funct);
				
		std::vector<real> dbls;	 //copy of every real value
		std::vector<char> ops;	 //copy of every operator in left-to-right order

		// Insert '+' signs before negative values.
		for (int i = 1; i < funct.length(); ++i)
		{
			if (funct[i] == '-' && std::isdigit(funct[i-1]))
			{
				funct.insert(funct.begin() + i, '+');
				i++;
			}
		}

		// Add real numbers to array.
		int idx = -1;
		while ((idx = findNearestReal(funct)) >= 0)
		{
			dbls.push_back(ParseInputNumber(getNearestReal(funct)));
			funct = removeNextReal(funct);
		}
		if (dbls.size() == 1)
			return dbls[0];

		// Collect leftover operators for arithmetic evaluation
		for (int i = 0; i < funct.length(); ++i) 
		{
			if ((isOperatorNotParen(funct[i]) && funct[i] != '-') ||
				(funct[i] == '-' &&
					((i == 0 && std::isdigit(funct[i + 1])) ||
						(i > 0 && funct[i - 1] != ',' && 
							funct[i - 1] != ')' && 
							!std::isdigit(funct[i - 1]))))) 
			{
				ops.push_back(funct[i]);
				funct.insert(funct.begin() + i, ',');
				funct.erase(i + 1, 1);
			}
		}
		funct.clear();

		if (ops.size() == 0) { return ParseInputNumber(funct); }//if operator vector is empty, there are no operators, only a value to return
		real answer = 0;
		int pass = 0;

		//evaluate in PEMDMAS order, 6 passes through the vector, eliminating evaluated terms and operators as you go
		while (dbls.size() > 1) {

			for (int i = 0; dbls.size() > 1 && i < ops.size(); ++i) {
				//note: format = 'a+b'

				real a = dbls[i];
				real b = dbls[i + 1];
				if (i >= 0 && ops[i] == '+' && pass == 4) {
					answer = a + b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '-' && pass == 5) {
					answer = a - b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '*' && pass == 1) {
					answer = a * b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '/' && pass == 2) {
					if (b == 0) {
						return NAN;//division by 0 is impossible
					}
					else {
						answer = a / b;
						dbls[i + 1] = answer;
						dbls.erase(dbls.begin() + i);
						ops.erase(ops.begin() + i);
						--i;
					}
				}
				if (i >= 0 && ops[i] == '^' && pass == 0) {
					answer = pow(a, b);
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
				if (i >= 0 && ops[i] == '%' && pass == 3) {
					answer = (int)a % (int)b;
					dbls[i + 1] = answer;
					dbls.erase(dbls.begin() + i);
					ops.erase(ops.begin() + i);
					--i;
				}
			}
			++pass;
			pass %= 6;//keep going around until only one dbl remains
		}
		return answer;
	}

	std::vector<real> ParseNInputNumbers(std::string str) {
		
		// Format input string to prevent issues during tokenization.
		str = removeSpaces(str);
		//str = replaceConstants(str);	//replace a constant like "PI" with value for pi.

		// First, make sure there are any numbers in this string to tokenize.
		std::vector<real> vals;
		bool shouldTokenize = false;
		for (int i = 0; i < str.length(); ++i)
		{
			if (std::isdigit(str[i]))
			{
				shouldTokenize = true;
				break;
			}
		}

		if (!shouldTokenize)
			return vals;

		// Tokenize string.
		std::vector<std::string> tokens = tokenize(str, ",");
		for (int i = 0; i < tokens.size(); ++i)
		{
			vals.push_back(ParseInputNumber(tokens[i]));
		}

		return vals;
	}


	complex ParseComplexNumber(std::string in) {
		//if the input string is only real-valued and not complex, evaluate it as a real number
		if (in.find('i') == std::string::npos) {
			return complex(ParseInputNumber(in), 0);
		}

		//else if it is complex...
		while (std::isdigit(in[0]) == false && in[0] != '-') { in = in.substr(1); }
		if (in.find('i') != std::string::npos) {
			if (isOperator(in[in.find('i') - 1])) //case: '1-i' is the number
			{
				in.insert(in.find('i'), "1");
			}
			if (in.find("+") == false && in.find("-") == false) {//case: '3i' is the number
				in = replaceChar(in, 'i', ' ');
				in.append(")");
				return complex(0, ParseInputNumber(in));
			}
		}

		for (int i = 1; i < in.size(); ++i) 
		{
			if (in[i] == '-' && std::isdigit(in[i - 1])) 			
				in.insert(in.begin() + i, ',');			
		}

		in = replaceChar(in, '+', ',');
		in = replaceChar(in, 'i', ' ');
		in.append(")");
		std::vector<real> temp = ParseNInputNumbers(in);
		if (temp.size() == 1) 
			return complex(temp[0], 0);
		return complex(temp[0], temp[1]);
	}
}