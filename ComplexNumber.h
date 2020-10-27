/**
*   Generic complex number class, written as a wrapper around the std::complex class.
*/

#pragma once
#include "Types.h"
#include "Vector.h"
#include <complex>
#include <tgmath.h>

namespace MathParser
{
	class ComplexNumber
	{
	private:
		// Wrapped std::complex object.
		std::complex<real> z;

	public:
		ComplexNumber(){}
		ComplexNumber(real Re) : z(Re, 0) {}
		ComplexNumber(real Re, real Im) : z(Re, Im) {}
		ComplexNumber(std::complex<real> z1) : z(z1) {}
		ComplexNumber(Vector v);

		real operator [](int i) const { return z._Val[i]; }
		real& operator [](int i) { return z._Val[0]; }
		bool operator==(const ComplexNumber& rhs){return z == rhs.z;}
		bool operator!=(const ComplexNumber& rhs){return !(z == rhs.z);}
		ComplexNumber& operator*=(ComplexNumber rhs);
		ComplexNumber& operator*=(real x);
		ComplexNumber& operator/=(ComplexNumber rhs);
		ComplexNumber& operator/=(real x);
		ComplexNumber& operator+=(ComplexNumber rhs);
		ComplexNumber& operator+=(real x);
		ComplexNumber& operator-=(ComplexNumber rhs);
		ComplexNumber& operator-=(real x);
		friend ComplexNumber operator+(ComplexNumber lhs, ComplexNumber rhs);
		friend ComplexNumber operator+(ComplexNumber lhs, real x);
		friend ComplexNumber operator+(real x, ComplexNumber rhs);
		friend ComplexNumber operator-(ComplexNumber lhs, ComplexNumber rhs);
		friend ComplexNumber operator-(ComplexNumber lhs, real x);
		friend ComplexNumber operator-(real x, ComplexNumber rhs);
		friend ComplexNumber operator*(ComplexNumber lhs, ComplexNumber rhs);
		friend ComplexNumber operator*(ComplexNumber lhs, real x);
		friend ComplexNumber operator*(real x, ComplexNumber rhs);
		friend ComplexNumber operator/(ComplexNumber lhs, ComplexNumber rhs);
		friend ComplexNumber operator/(ComplexNumber lhs, real x);
		friend ComplexNumber operator/(real x, ComplexNumber rhs);

		real im() { return z.imag(); }
		real re() { return z.real(); }

	private:
		// Helper functions to reduce reusing code for operator overloads.
		ComplexNumber add(ComplexNumber a, ComplexNumber b);
		ComplexNumber add(ComplexNumber a, real  b);
		ComplexNumber add(real b, ComplexNumber a);
		ComplexNumber subtract(ComplexNumber a, ComplexNumber b);
		ComplexNumber subtract(ComplexNumber a, real  b);
		ComplexNumber subtract(real b, ComplexNumber a);
		ComplexNumber multiply(ComplexNumber a, ComplexNumber b);
		ComplexNumber multiply(ComplexNumber a, real b);
		ComplexNumber multiply(real b, ComplexNumber a);
		ComplexNumber divide(ComplexNumber p, ComplexNumber q);
		ComplexNumber divide(ComplexNumber p2, real q);
		ComplexNumber divide(real q, ComplexNumber p2);
	};
}

namespace std
{
	// Wrapper for standard 'pow()' function to use ComplexNumber class.
	MathParser::ComplexNumber pow(MathParser::ComplexNumber z, MathParser::ComplexNumber b);
	template <class T> MathParser::ComplexNumber pow(MathParser::ComplexNumber z, T b);
}