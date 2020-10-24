#include "Function.h"
#include "StringUtils.h"
#include "Parsing.h"
#include "MathLib.h"

namespace MathParser
{
	Function::Function(std::string str) {		//input of format f(x,y,w,z...) = x*3 + exp(-3*y*z) - ...
		if (str.find("=") != std::string::npos) {	//case:  string is 'f(x) = 2*x + exp(-y)'
			while (str[0] != '(') { str = str.substr(1); }//clip off everything leading up to the variable definition
			str = str.substr(1);
			for (int i = 0; i < str.find("=") + 1; ++i) {
				if (std::isalpha(str[i])) 
					variables.push_back(str[i]);
			}
			function = str.substr(str.find("=") + 1);
		}

		else { function = str; }				//case:  string has no variables, i.e. '23-12+exp(3)'
	}

	Function::Function(std::string vars, std::string funct) {
		for (int i = 0; i < vars.size(); ++i) {
			if (std::isalpha(vars[i])) 
				variables.push_back(vars[i]);
		}
		function = funct;
	}

	Function::Function(std::vector<char> vars, std::string functs) {
		function = functs;
		variables = vars;
	}

	Function::Function(std::vector<std::string> vec)
	{
		function = vec[1];
		variables = varsToChars(vec[0]);
	}

	std::vector<char> Function::varsToChars(std::string vars) {
		std::vector<char> variables;
		for (int i = 0; i < vars.size(); ++i)
			if (std::isalpha(vars[i]) && vars[i] != ',') 
				variables.push_back(vars[i]);
		return variables;
	}

	std::string Function::to_string(int precision)
	{
		std::string str = "";
		if (!variables.size() && !function.length())
			return str;

		str = "f(";
		for (int i = 0; i < variables.size(); ++i)
		{
			std::string temp = " ";
			temp[0] = variables[i];
			str.append(temp);

			if (i < variables.size() - 1)
				str.append(",");
		}
		str.append(") = ");
		str.append(function);

		return str;
	}

	real Function::evalFunction(std::string str) 
	{
		// User-defined function names should have 2 or more letters in a row, else
		// they are likely to be variables.
		if (str.length() < 2)
			return 0; 

		else if(str.find("abs(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("abs(") + 1));
			return std::abs(tmp);
		}
		else if(str.find("acos(") != std::string::npos && str.find("acosh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("acos(") + 1));
			return acos(tmp);
		}
		else if(str.find("acosh(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("acosh(") + 1));
			return acosh(tmp);
		}
		else if(str.find("asin(") != std::string::npos && str.find("asinh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("asin(") + 1));
			return asin(tmp);
		}
		else if(str.find("asinh(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("asinh(") + 1));
			return asinh(tmp);
		}
		else if(str.find("atan(") != std::string::npos && str.find("atan2(") == std::string::npos && str.find("atanh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("atan(") + 1));
			return atan(tmp);
		}
		else if(str.find("atan2(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return atan2(tmp[1], tmp[0]);
		}
		else if(str.find("atanh(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("atanh(") + 1));
			return atanh(tmp);
		}
		else if(
			str.find("BesselFunctionFirstKind(") != std::string::npos
			&& str.find("sphericalBesselFunctionFirstKind") == std::string::npos
			&& str.find("regularModifiedBesselFunctionFirstKind(") == std::string::npos
			) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return std::cyl_bessel_j(tmp[0], tmp[1]);
		}
		else if(str.find("BesselFunctionSecondKind(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return BesselFunction2ndKind(tmp[0], tmp[1]);
		}
		else if(str.find("beta(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return std::beta(tmp[0], tmp[1]);
		}
		else if(str.find("cbrt(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("cbrt(") + 1));
			return std::cbrt(tmp);
		}
		else if(str.find("combination(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("combination(") + 1));
			if(tmp.size() == 2) 
				return combination(tmp[0], tmp[1]);
			else 
				return 0;
		}
		else if(str.find("ceil(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("ceil(") + 1));
			return ceil(tmp);
		}
		else if(str.find("cos(") != std::string::npos && str.find("acos(") == std::string::npos && 
			str.find("cosh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("cos(") + 1));
			return cos(tmp);
		}
		else if(str.find("cosh(") != std::string::npos && str.find("acosh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("cosh(") + 1));
			return cosh(tmp);
		}
		else if(str.find("cot(") != std::string::npos && str.find("coth(") == std::string::npos && str.find("acot(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("cot(") + 1));
			return 1.0/tan(tmp);
		}
		else if(str.find("csc(") != std::string::npos && str.find("acsc(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("csc(") + 1));
			return 1.0/sin(tmp);
		}
		else if(str.find("digamma(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return digamma(tmp);
		}
		else if(str.find("doubleFactorial(") != std::string::npos)
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return doubleFactorial(tmp);
		}
		else if(str.find("erf(") != std::string::npos)
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return std::erf(tmp);
		}
		else if(str.find("erfc(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("erfc(") + 1));
			return std::erfc(tmp);
		}
		else if(str.find("exp(") != std::string::npos && str.find("exp2(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("exp(") + 1));
			return exp(tmp);
		}
		else if(str.find("exp2(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("exp2(") + 1));
			return exp2(tmp);
		}
		else if(str.find("factorial(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("factorial(") + 1));
			return factorial(tmp);
		}
		else if(str.find("fallingFactorial(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return fallingFactorial(tmp[0], tmp[1]);
		}
		else if(str.find("floor(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("floor(") + 1));
			return std::floor(tmp);
		}
		else if(str.find("gamma(") != std::string::npos &&
			str.find("digamma(") == std::string::npos)
		{
			real tmp = ParseInputNumber(str.substr(str.find("gamma(") + 1));
			return gamma(tmp);
		}
		else if(str.find("hypot(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return sqrt(pow(tmp[0], 2) + pow(tmp[1], 2));
		}
		else if(str.find("int(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("int(") + 1));
			return (int)(tmp);
		}
		else if(str.find("LambertW(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return LambertW(tmp);
		}
		else if(str.find("log(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			if(tmp.size() == 1) //standard log base e
				return log(tmp[0]);
			else if(tmp.size() == 2)//log w/different base
				return logarithm(tmp[0],tmp[1]);
			return 0;
		}
		else if(str.find("log2(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("log2(") + 1));
			return log2(tmp);
		}
		else if(str.find("log10(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("log10(") + 1));
			return log10(tmp);
		}
		else if(str.find("logGamma(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return std::lgamma(tmp);
		}
		else if(str.find("logit(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("logit(") + 1));
			return logit(tmp);
		}
		else if(str.find("max(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			if (tmp.size() == 2)
				return (tmp[0] > tmp[1]) ? tmp[0] : tmp[1];
			else
				return 0;
		}
		else if(str.find("min(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			if (tmp.size() == 2)
				return (tmp[0] < tmp[1]) ? tmp[0] : tmp[1];
			else 
				return 0;
		}
		else if(str.find("permutation(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			if(tmp.size() == 2) 
				return permutation(tmp[0], tmp[1]);
			else 
				return 0;
		}
		else if(str.find("pow(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return pow(tmp[0], tmp[1]);
		}
		else if(str.find("prime(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return prime(tmp);
		}
		else if(str.find("sumPrimes(") != std::string::npos) 
		{
			int tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return sumPrimes(tmp);
		}
		else if(str.find("productPrimes(") != std::string::npos) 
		{
			int tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return productPrimes(tmp);
		}
		else if(str.find("primeCountingFunction(") != std::string::npos) 
		{
			int tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return primeCountingFunction(tmp);
		}
		else if(str.find("rect(") != std::string::npos) {
			std::string s = str.substr(str.find("rect("));
			s = s.substr(str.find("(") + 1, str.find(")") - 1);
			real tmp = ParseInputNumber(s);
			return rect(tmp);
		}
		else if(str.find("regularModifiedBesselFunctionFirstKind(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return std::cyl_bessel_i(tmp[0], tmp[1]);
		}
		else if(str.find("risingFactorial(") != std::string::npos) 
		{
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return risingFactorial(tmp[0], tmp[1]);
		}
		else if(str.find("sec(") != std::string::npos 
			&& str.find("asec(") == std::string::npos 
			&& str.find("sech(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("sec(") + 1));
			return 1.0/cos(tmp);
		}
		else if (str.find("sech(") != std::string::npos)
		{
			real tmp = ParseInputNumber(str.substr(str.find("sech(") + 1));
			return 1.0/cosh(tmp);
		}
		else if (str.find("sigmoid(") != std::string::npos)
		{
			real tmp = ParseInputNumber(str.substr(str.find("sigmoid(") + 1));
			return sigmoid(tmp);
		}
		else if(str.find("sin(") != std::string::npos && str.find("asin(") == std::string::npos && str.find("sinh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("sin(") + 1));
			return sin(tmp);
		}
		else if(str.find("sinc(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return sinc(tmp);
		}
		else if(str.find("sinh(") != std::string::npos && str.find("asinh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("sinh(") + 1));
			return sinh(tmp);
		}
		else if(str.find("sgn(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("sgn(") + 1));
			return sgn(tmp);
		}
		else if(str.find("sphericalBesselFunctionFirstKind(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return std::sph_bessel(tmp[0], tmp[1]);
		}
		else if(str.find("sqrt(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("sqrt(") + 1));
			return sqrt(tmp);
		}
		else if(str.find("squarewave(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("squarewave(") + 1));
			return squarewave(tmp);
		}
		else if(str.find("step(") != std::string::npos) {//Heaviside step function
			std::string s = str.substr(str.find("step("));
			s = s.substr(str.find("(") + 1, str.find(")") - 1);
			real tmp = ParseInputNumber(s);
			return step(tmp);
		}
		else if(str.find("StirlingNumberFirstKind(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return StirlingNumber1stKind(tmp[0], tmp[1]);
		}
		else if(str.find("StirlingNumberSecondKind(") != std::string::npos) {
			std::vector<real> tmp = ParseNInputNumbers(str.substr(str.find("(") + 1));
			return StirlingNumber2ndKind(tmp[0], tmp[1]);
		}
		else if(str.find("tan(") != std::string::npos && str.find("tanh(") == std::string::npos && str.find("atan(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("tan(") + 1));
			return tan(tmp);
		}
		else if(str.find("tanh(") != std::string::npos && str.find("atanh(") == std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("tanh(") + 1));
			return tanh(tmp);
		}
		else if(str.find("tri(") != std::string::npos) {
			std::string s = str.substr(str.find("tri("));
			s = s.substr(str.find("(") + 1, str.find(")") - 1);
			real tmp = ParseInputNumber(s);
			return tri(tmp);
		}
		else if(str.find("trunc(") != std::string::npos) {
			real tmp = ParseInputNumber(str.substr(str.find("acos(") + 1));
			return std::trunc(tmp);
		}
		else if(str.find("WrightOmegaFunction(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("WrightOmegaFunction(") + 1));
			return WrightOmegaFunction(tmp);
		}
		else if(str.find("zeta(") != std::string::npos) 
		{
			real tmp = ParseInputNumber(str.substr(str.find("(") + 1));
			return std::riemann_zeta(tmp);
		}

		return 0;	//else, if none of the other functions match, return NULL
	}

	real Function::evaluate(std::vector<real> vals, int precision) {
		/*	NOTE:
		Suitable arithmetic operators are +,-,*,/,^.  Suitable functions are found in evalFunction method.
		This function will properly evaluate a function of the form f(x,y) = x^2 without error.		*/

		/*=================
		|DECLARE VARIABLES|
		=================*/
		//make sure when using to_string precision or other functions which truncate values, we do as little truncation as possible
		int prcs = precision;

		std::string funct = function;//local copy of function string, necessary as we will erase/replace parts of the string as we go

		funct = removeSpaces(funct);//remove all spaces in function string
		//funct = replaceConstants(funct);//replace all system constants in string
		//========================================

		/*================================
		|REPLACE VARIABLES WITH CONSTANTS|
		================================*/
		//replace all variables with their real value from real* vals
		if (vals.size() > 0) {
			for (int i = 0; i < variables.size(); ++i) {
				std::string val = to_stringPrecision(vals[i], precision);//real values for variable made into a string

			   //check string for spots to replace char, then perform replace function
				for (int j = 0; j < funct.size(); ++j) {
					if (j == 0 && funct[j] == variables[i] && std::isalpha(funct[j + 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
						j = 0;
					}
					if (j > 0 && j < funct.size() - 1 && funct[j] == variables[i] && std::isalpha(funct[j + 1]) == false && std::isalpha(funct[j - 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
					}
					if (j == funct.size() - 1 && funct[j] == variables[i] && std::isalpha(funct[j - 1]) == false) {
						funct.erase(j, 1);
						funct.insert(j, val);
					}
				}
			}
		}//================================================================

		 //remove any '-' signs that might have resulted accidentally
		for (int i = 0; i < funct.size() - 1; ++i) {
			if (funct[i] == '-' && funct[i + 1] == '-') {
				funct[i] = ' ';
				funct[i + 1] = '+';
				++i;
			}
		}
		funct = removeSpaces(funct);

		/*================
		|HANDLE FUNCTIONS|
		================*/
		bool searchForFuncts = hasFunction(funct);
		while (searchForFuncts) 
		{
			bool foundOne = false;
			for (int i = 0; i < funct.size() - 1; ++i) {
			topOfTheFunctsLoop:
				if (std::isalpha(funct[i]) && std::isalpha(funct[i + 1])) {		//find functions, which have more than one letter in a row
					std::string tempstr = funct.substr(i);
					tempstr = tempstr.substr(0, findMatchingParen(tempstr, tempstr.find("(")) + 1);		//save function (i.e. 'exp(12.2)') as a substr to evaluate
					if (hasFunction(tempstr.substr(tempstr.find("(")))) {//keep searching for the innermost paren without a function inside
						i += tempstr.find("(") - 1;//skip past the function (i.e. i+=4 if outer function is 'cos(' for example)
						goto topOfTheFunctsLoop;
					}
					std::string toReplace = tempstr;

					// Make sure there's nothing within the parenthesis that needs to be 
					// simplified (i.e. 'sin(12-2)').
					int posit = tempstr.find('(');
					if (hasOperator(tempstr.substr(posit + 1, findMatchingParen(tempstr, posit) + 1)) 
						&& countReals(tempstr) > 1) {
						int opers = 0;
						for (int k = 0; k < tempstr.length() - 2; ++k) {
							if (isOperator(tempstr[k])) {
								if (tempstr[k] == '-') {
									if (k > 0) {
										if (std::isdigit(tempstr[k - 1]) && 
											std::isdigit(tempstr[k + 1])) {
											++opers;
										}
									}
								}
								else { ++opers; }
							}
						}

						if (opers > 0) {
							size_t a = tempstr.find("(") + 1;
							size_t b = tempstr.find(")") - 1;
							std::string chk = tempstr.substr(a, b);
							chk = replaceChar(chk, ')', ' ');
							chk = replaceChar(chk, '(', ' ');
							chk = removeSpaces(chk);
							Function f2(chk);
							real p = f2.evaluate(0.0);
							std::string replc = to_stringPrecision(p, precision);
							tempstr = replaceString(tempstr, chk, replc);
						}
					}

					//now evaluate the function, get the real value, and replace the whole thing in the original string (i.e. 'sin(0)' --> '0')
					real tmpdbl = evalFunction(tempstr);
					funct = replaceString(funct, toReplace, to_stringPrecision(tmpdbl, precision));
					foundOne = true;
				}
			}
			if (foundOne == false) { searchForFuncts = false; }
		}//====================================================
		if (countReals(funct) == 1) {
			precision = prcs;
			return ParseInputNumber(getNearestReal(funct));
		}
		funct = processInnerTerms(funct); //handle parens
		real answer = calculateArithmetic(funct);//calculate arithmetic of the last few values
		precision = prcs; //finally, restore precision and return a value
		return answer; //note: if unsuccessful, will return 0
	}//END EVALUATE FUNCTION================================================

	real Function::evaluate(real x, int precision) {
		std::vector<real> vals;
		vals.push_back(x);
		real n = evaluate(vals,precision);
		return n;
	}

}