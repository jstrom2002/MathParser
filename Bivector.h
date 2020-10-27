#pragma once
#include "Types.h"
#include "Vector.h"

namespace MathParser
{
	class Bivector
	{
	public:
		Bivector(){}
		Bivector(Vector a_, Vector b_) : a(a_), b(b_) {}

		// Overloaded operators.
		Vector operator [](int i) const;
		Vector& operator [](int i);
		bool operator==(const Bivector& rhs) { return (a==rhs.a) && (b==rhs.b); }
		bool operator!=(const Bivector& rhs) { return !(*this==rhs); }
		Bivector& operator*=(Bivector rhs);
		Bivector& operator*=(real x);
		Bivector& operator/=(Bivector rhs);
		Bivector& operator/=(real x);
		Bivector& operator+=(Bivector rhs);
		Bivector& operator+=(real x);
		Bivector& operator-=(Bivector rhs);
		Bivector& operator-=(real x);
		friend Bivector operator+(Bivector& lhs, Bivector& rhs);
		friend Bivector operator+(Bivector& lhs, real x);
		friend Bivector operator+(real x, Bivector& rhs);
		friend Bivector operator-(Bivector& lhs, Bivector& rhs);
		friend Bivector operator-(Bivector& lhs, real x);
		friend Bivector operator-(real x, Bivector& rhs);
		friend Bivector operator*(Bivector lhs, Bivector& rhs);
		friend Bivector operator*(Bivector& lhs, real x);
		friend Bivector operator*(real x, Bivector& rhs);
		friend Bivector operator/(Bivector& lhs, Bivector& rhs);
		friend Bivector operator/(Bivector& lhs, real x);
		friend Bivector operator/(real x, Bivector& rhs);

	private:
		// Wrapped 'Vector' objects representing the results of exterior/wedge product.
		Vector a, b;

		// Helper functions to reduce reusing code for operator overloads.
		Bivector add(Bivector& a, Bivector& b);
		Bivector add(Bivector& a, real  b);
		Bivector add(real b, Bivector& a);
		Bivector subtract(Bivector& a, Bivector& b);
		Bivector subtract(Bivector& a, real b);
		Bivector subtract(real b, Bivector& a);
		Bivector multiply(Bivector& a, Bivector& b);
		Bivector multiply(Bivector& a, real b);
		Bivector multiply(real b, Bivector& a);
		Bivector multiply(Bivector& a, Vector& b);
		Bivector multiply(Vector b, Bivector& a);
		Bivector divide(Bivector& p, Bivector& q);
		Bivector divide(Bivector& p2, real q);
		Bivector divide(real q, Bivector& p2);
	};
}