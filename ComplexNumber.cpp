#include "ComplexNumber.h"

namespace std
{
	MathParser::ComplexNumber pow(MathParser::ComplexNumber z, MathParser::ComplexNumber b)
	{
		return std::pow(std::complex(z.re(), z.im()), std::complex(b.re(), b.im()));
	}

	MathParser::ComplexNumber pow(MathParser::ComplexNumber z, MathParser::real b)
	{
		return std::pow(std::complex(z.re(), z.im()), b);
	}
}

namespace MathParser
{
	ComplexNumber ComplexNumber::add(ComplexNumber a, ComplexNumber b){	return ComplexNumber(a.z + b.z);}
	ComplexNumber ComplexNumber::add(ComplexNumber a, real  b){return ComplexNumber(a.z + b);}
	ComplexNumber ComplexNumber::add(real b, ComplexNumber a){return ComplexNumber(a.z + b);}
	ComplexNumber ComplexNumber::subtract(ComplexNumber a, ComplexNumber b){return ComplexNumber(a.z - b.z);}
	ComplexNumber ComplexNumber::subtract(ComplexNumber a, real  b){return ComplexNumber(a.z - b);}
	ComplexNumber ComplexNumber::subtract(real b, ComplexNumber a){return ComplexNumber(a.z - b);}
	ComplexNumber ComplexNumber::multiply(ComplexNumber a, ComplexNumber b){return ComplexNumber(a.z * b.z);}
	ComplexNumber ComplexNumber::multiply(ComplexNumber a, real b){	return ComplexNumber(a.z * b);}
	ComplexNumber ComplexNumber::multiply(real b, ComplexNumber a){return ComplexNumber(a.z * b);}
	ComplexNumber ComplexNumber::divide(ComplexNumber a, ComplexNumber b){return ComplexNumber(a.z / b.z);}
	ComplexNumber ComplexNumber::divide(ComplexNumber a, real b){return ComplexNumber(a.z / b);}
	ComplexNumber ComplexNumber::divide(real b, ComplexNumber a){return ComplexNumber(a.z / b);}
	ComplexNumber& ComplexNumber::operator*=(ComplexNumber rhs) { *this = *this * rhs; return *this; }
	ComplexNumber& ComplexNumber::operator*=(real x) { *this = *this * x; return *this; }
	ComplexNumber& ComplexNumber::operator/=(ComplexNumber rhs) { *this = *this / rhs; return *this; }
	ComplexNumber& ComplexNumber::operator/=(real x) { *this = *this / x; return *this; }
	ComplexNumber& ComplexNumber::operator+=(ComplexNumber rhs) { *this = *this + rhs; return *this; }
	ComplexNumber& ComplexNumber::operator+=(real x) { *this = *this + x; return *this; }
	ComplexNumber& ComplexNumber::operator-=(ComplexNumber rhs) { *this = *this - rhs; return *this; }
	ComplexNumber& ComplexNumber::operator-=(real x) { *this = *this - x; return *this; }
	ComplexNumber operator+(ComplexNumber lhs, ComplexNumber rhs) { return lhs.add(lhs, rhs); }
	ComplexNumber operator+(ComplexNumber lhs, real x) { return lhs.add(lhs, x); }
	ComplexNumber operator+(real x, ComplexNumber rhs) { return rhs.add(x, rhs); }
	ComplexNumber operator-(ComplexNumber lhs, ComplexNumber rhs) { return lhs.subtract(lhs, rhs); }
	ComplexNumber operator-(ComplexNumber lhs, real x) { return lhs.subtract(lhs, x); }
	ComplexNumber operator-(real x, ComplexNumber rhs) { return rhs.subtract(x, rhs); }
	ComplexNumber operator*(ComplexNumber lhs, ComplexNumber rhs) { return lhs.multiply(lhs, rhs); }
	ComplexNumber operator*(ComplexNumber lhs, real x) { return lhs.multiply(lhs, x); }
	ComplexNumber operator*(real x, ComplexNumber rhs) { return rhs.multiply(x, rhs); }
	ComplexNumber operator/(ComplexNumber lhs, ComplexNumber rhs) { return lhs.divide(lhs, rhs); }
	ComplexNumber operator/(ComplexNumber lhs, real x) { return lhs.divide(lhs, x); }
	ComplexNumber operator/(real x, ComplexNumber rhs) { return rhs.divide(x, rhs); }
}