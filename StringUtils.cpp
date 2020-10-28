#pragma once
#include "StringUtils.h"
#include "ComplexNumber.h"
#include <sstream>
#include <iomanip>
#include <algorithm>

namespace MathParser
{
	bool isInteger(std::string n) {
		if (n.find('.') == std::string::npos) { return true; }
		return false;
	}

	std::string getMantissa(std::string n) {
		if (isInteger(n)) 
			return "";
		return n.substr(n.find('.') + 1, n.length() - n.find('.') - 1);
	}

	std::string getInteger(std::string n) {
		if (isInteger(n)) 
			return n;
		return n.substr(0, n.find('.'));
	}

	bool isOperator(char c) {//valid operators are +,-,*,/,^,(,)
		if (
			c == '+' ||
			c == '-' ||
			c == '*' ||
			c == '/' ||
			c == '^' ||
			c == '%' ||
			c == '(' ||
			c == ')'
			) {
			return true;
		}
		else { return false; }
	}

	bool isOperatorNotParen(char c) {//operators are +,-,*,/,^
		if (
			c == '+' ||
			c == '-' ||
			c == '*' ||
			c == '/' ||
			c == '^' ||
			c == '%'
			) {
			return true;
		}
		else { return false; }
	}

	char randomChar() {
		srand(time(0) + clock());
		char c = (rand() % 255) + 33;
		while (c == 127) { c = (rand() % 255) + 33; }
		return c;
	}

	std::string replaceChar(std::string str, char replace, char replacer) {
		std::replace(str.begin(), str.end(), replace, replacer);
		return str;
	}

	std::string replaceString(std::string str, std::string replace, std::string replacer) {
		for (int i = 0; i <= str.length(); ++i) {//run replace loop ONCE only
			std::string sub = str.substr(i, str.length());
			int pos = sub.find(replace) + i;
			if (sub.find(replace) == std::string::npos || pos < 0 || sub.find(replace) > str.length() - 1) {
				//break if done replacing
				i = str.length() + 1;
				pos = std::string::npos;//set position to NULL for string type
			}
			if (pos >= 0) {
				str.erase(pos, replace.length());
				str.insert(pos, replacer);
				i = pos + replacer.length()-1;
			}
		}
		return str;
	}

	std::string removeSpaces(std::string str) 
	{
		std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
		str.erase(end_pos, str.end());
		return str;
	}

	std::string trim(std::string str) //removes spaces in front/back of string
	{
		while (str[str.length() - 1] == ' ') { str.pop_back(); }
		while (str[0] == ' ') { str = str.substr(1); }
		return str;
	}

	std::string reverseString(std::string str) 
	{
		std::string s = str;
		std::reverse(s.begin(), s.end());
		return s;
	}

	std::string to_string_precision(real val, unsigned int precisn) 
	{
		std::ostringstream out;
		int p = precisn;
		if (val == floor(val)) { p = 0; }
		out << std::fixed;
		out << std::setprecision(p) << val;
		out << std::scientific;
		return out.str();
	}

	std::string to_string_precision(ComplexNumber z, unsigned int precision)
	{
		std::string answer;
		std::ostringstream strs;
		int prc = precision;
		if (z.re() == std::floor(z.re()))
			prc = 0;
		strs << std::fixed << std::setprecision(prc) << z.re();
		answer.append(strs.str());
		prc = precision;
		if (z.im() == floor(z.im()))
			prc = 0;
		if (z.im() != 0) {
			if (z.im() < 0) {
				answer.append("-");
			}
			else {
				answer.append("+");
			}
			if (std::abs(z.im()) != 1) {
				std::ostringstream strs2;
				strs2 << std::fixed << std::setprecision(prc) << std::abs(z.im());
				answer.append(strs2.str());
			}
			answer.append("i");
		}
		return answer;
	}

	std::string toHex(int val) {
		std::stringstream stream;
		stream << "0x";
		if (val < 16) {
			stream << "0";
		}
		stream << std::hex << std::uppercase << val;
		return stream.str();
	}

	std::string to_hex(std::string input) {
		int val = std::stod(input);
		return toHex(val);
	}

	std::string matchParens(std::string str) {
		int frontParens = 0;
		int backParens = 0;
		for (int i = 0; i <= str.length(); ++i) {
			if (str[i] == '(') { ++frontParens; }
			if (str[i] == ')') { ++backParens; }
		}
		while (frontParens > backParens) {
			str.append(")");
			++backParens;
		}
		return str;
	}

	size_t findMatchingParen(std::string str, int pos) {
		/*given the position 'pos' of a front paren at a position in string 'str',
		this function returns the position of the matching closing paren (or 
		std::string::npos if it does not exist). */
		int backParens = 0;
		if (str[pos] == '(') {//case: str[pos] is an open bracket
			for (int i = pos + 1; i <= str.length(); ++i) {
				if (str[i] == '(') { ++backParens; }
				if (str[i] == ')') { --backParens; }
				if (backParens < 0) { return i; }
			}
		}
		if (str[pos] == ')') {//case: str[pos] is an closing bracket
			for (int i = pos - 1; i >= 0; --i) {
				if (str[i] == '(') { --backParens; }
				if (str[i] == ')') { ++backParens; }
				if (backParens < 0) { return i; }
			}
		}
		return std::string::npos;//else, return null position
	}

	std::string removeUnmatchedParens(std::string funct) {
		for (int i = 0; i < funct.length(); ++i) {
			if (funct[i] == '(' && findMatchingParen(funct, i) == std::string::npos) {
				funct.erase(i, 1);
			}
		}
		return funct;
	}

	std::string eraseParenPair(std::string str, size_t pos) {
		str.erase(findMatchingParen(str, pos), 1);
		str.erase(pos, 1);
		return str;
	}

	std::vector<std::string> tokenize(std::string toTokenize, std::string token) 
	{
		std::vector<std::string> result;
		char* tk = strtok((char*)toTokenize.c_str(), token.c_str());

		// Keep printing tokens while one of delimiters are present. 
		while (tk != NULL)
		{
			result.push_back(tk);
			tk = strtok(NULL, token.c_str());
		}

		return result;
	}

	std::string getDirectory(std::string str) 
	{
		if (str.find("/") != std::string::npos) 
			str = str.substr(0, str.rfind("/") + 1);
		if (str.find("\\") != std::string::npos)	
			str = str.substr(0, str.rfind("\\") + 1);
		return "";
	}

	std::string getExtension(std::string str) 
	{// File extensions will be converted to lowercase letters.
		std::string str2 = str;
		std::transform(str2.begin(), str2.end(), str2.begin(), std::tolower);
		if (str.find(".") == std::string::npos)
			return "";
		str = str.substr(str.rfind("."));
		return str;
	}

	std::string getFilename(std::string str) 
	{
		while (str.rfind("/") != std::string::npos && str.length() > 0)
			str = str.substr(str.rfind("/") + 1);		
		while (str.rfind("\\") != std::string::npos && str.length() > 0) 
			str = str.substr(str.rfind("\\") + 1);		
		return str;
	}
}