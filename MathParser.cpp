#include "MathParser.h"
#include "StringUtils.h"
#include "Polynomial.h"
#include "Function.h"
#include "Parsing.h"
#include <string>
#include <algorithm>

namespace MathParser
{
	std::string MathParser::printSavedFunctions()
	{
		std::string str = "";
		if (!savedFunctions.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < savedFunctions.size(); ++i)
		{
			std::string funcName = savedFunctions[i].to_string();
			funcName.insert(funcName.begin() + 1, '0' + i);
			str.append(funcName + newline);
		}

		return str;
	}

	std::string MathParser::printSavedPolynomials()
	{
		std::string str = "";
		if (!savedPolynomials.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < savedPolynomials.size(); ++i)
		{
			std::string funcName = "p(x) = ";
			funcName.insert(funcName.begin() + 1, '0' + i);
			funcName.append(savedPolynomials[i].to_string());
			str.append(funcName + newline);
		}

		return str;
	}

	std::string MathParser::parse(std::string str)
	{
		if (!str.length() || str == " " || str == "\n") 
			return ""; 

		// Check for functions, polynomials, etc.
		if (str.find("=") != std::string::npos)
		{
			if (str.find("p(") != std::string::npos)			
				savedPolynomials.push_back(Polynomial(ParsePolynomial(str)));			
			else			
				savedFunctions.push_back(Function(ParseFunction(str)));			
			return "";
		}

		// If no functions, then parse arithmetic string.
		if (hasOperator(str)) 
		{
			// Remove all spaces in input to format input properly.
			str = removeSpaces(str);

			Function fp("", str);
			return std::to_string(fp.evaluate(0, this->precision));
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

		//SEPARATE VALUES
		//===============
		//cut up input string into smaller strings and place them in the "vals" string vector
		std::vector<std::string> vals;
		for (int i = 0; i < input.length(); ++i) {
			if (i > 0 && input[i] == '+') {
				std::string s = input.substr(0, i);
				if (s[0] == '+') { s[0] = ' '; }
				removeSpaces(s);
				while (s[0] == ' ') { s = s.substr(1); }//makes sure to remove blank space at the beginning of the string
				vals.push_back(s);
				input.erase(0, i);
				i = 0;
			}
		}
		//deal with the last of the input string
		if (input[0] == '+') { input[0] = ' '; }
		removeSpaces(input);
		while (input[0] == ' ') { input = input.substr(1); }//makes sure to remove blank space at the beginning of the string
		vals.push_back(input);


		//GET MATCHING EXPONENT-COEFFICIENT PAIRS
		//=======================================
		//now cut the strings in the vector into two parts:  the coefficient and the exponent
		std::vector<std::string> indices;
		std::vector<std::string> coefs;
		for (int i = 0; i < vals.size(); ++i) {
			if (vals[i].find("^") == std::string::npos) {
				if (vals[i].find("x") != std::string::npos) //case: no "^" but an x.  i.e. '23.23x'
				{
					indices.push_back("1");
					std::string temp = vals[i];
					bool onlyX = true;
					for (int j = 0; j < temp.length(); ++j) {//check to see if it's just 'x' or if it has a coefficient (i.e. '2.34x')
						if (std::isdigit(temp[j])) 
							onlyX = false;
					}
					if (onlyX) //put a '1' at the beginning if there isn't one
						temp.insert(temp.begin(), '1');				
					std::replace(temp.begin(), temp.end(), 'x', ' ');
					removeSpaces(temp);
					coefs.push_back(temp);
				}
				if (vals[i].find("x") == std::string::npos) {//case: no o'^', no 'x'.  A single value alone (i.e. '12.32x^0')
					indices.push_back("0");
					coefs.push_back(vals[i]);
				}
			}
			else {
				for (int j = 0; j < vals[i].size(); ++j) {//case" a value with an 'x' and a '^' which needs separating in to coefficient and exponent
					if (j > 0 && vals[i][j] == '^') {
						std::string s = vals[i].substr(j + 1);
						std::string c = vals[i].substr(0, j - 1);
						if (s[0] == '^') 
							s[0] = ' ';
						removeSpaces(s);
						
						//makes sure to remove blank space at the beginning of the string
						while (s[0] == ' ') 
							s = s.substr(1);
						indices.push_back(s);

						//check now to see if there actually is a coefficient or if we've got a case like 'x^4'
						bool isCoefficient = false;
						for (int q = 0; q < c.size(); ++q) 
						{ 
							if (std::isdigit(c[q])) 
							{ 
								isCoefficient = true;	
								q = c.size(); 
							} 
						}

						//if there are no numbers, put a "1" string in there
						if (!isCoefficient) 
							coefs.push_back("1");
						else //else, treat normally
							coefs.push_back(c);
					}
				}
			}
		}

		//convert all strings to real/int values
		//======================================
		int largestExponent = ParseInputNumber(indices[0]);
		++largestExponent;
		std::vector<real> v;
		v.resize(largestExponent + 1);
		for (int i = 0; i < largestExponent + 1; ++i) {
			v[i] = 0;
		}
		for (int i = 0; i < indices.size(); ++i) {
			v[(int)ParseInputNumber(indices[i])] = ParseInputNumber(coefs[i]);
		}
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
		if (position.size() == 0) { answer.push_back(Function(ParseFunction(str))); return answer; }//if there is only one single value in the string
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
}