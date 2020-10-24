#include "complex.h"
#include "Types.h"
#include <sstream>
#include <iomanip>

namespace MathParser
{
	std::string toString(complex z, int precision) {
		std::string answer;
		std::ostringstream strs;
		int prc = precision;
		if (z.real() == std::floor(z.real())) 
			prc = 0;
		strs << std::fixed << std::setprecision(prc) << z.real();
		answer.append(strs.str());
		prc = precision;
		if (z.imag() == floor(z.imag()))
			prc = 0;
		if (z.imag() != 0) {
			if (z.imag() < 0) {
				answer.append("-");
			}
			else {
				answer.append("+");
			}
			if (std::abs(z.imag()) != 1) {
				std::ostringstream strs2;
				strs2 << std::fixed << std::setprecision(prc) << std::abs(z.imag());
				answer.append(strs2.str());
			}
			answer.append("i");
		}
		return answer;
	}

	std::string toStringBothParts(complex z, int precision) {
		std::string str = toString(z, precision);
		if (z.imag() == 0) {
			str.append("+0i");
		}
		return str;
	}
}