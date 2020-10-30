/**
*   Generic vector class, written as a wrapper around the std::vector class.
*/

#pragma once
#include "Types.h"
#include "BaseObject.h"
#include <vector>
#include <iterator>

namespace MathParser
{
	class ComplexNumber;
	class Matrix;
	class Polynomial;

	class Vector : public BaseObject
	{
	public:
		Vector(){}
		Vector(int i) { v.resize(i, 0); }
		Vector(std::vector<real> v);
		Vector(Polynomial p);
		Vector(ComplexNumber p);
		
		// Overloads of std namespace functions for std::vector class,
		// written for brevity (ie users can use 'Vector::begin()' instead
		// of 'Vector::v::begin()').
		std::vector<real>::iterator begin() { return v.begin(); }
		void clear() { v.clear(); }
		std::vector<real>::iterator end() { return v.end(); }
		void erase(std::vector<real>::iterator beg){v.erase(beg);}
		void erase(std::vector<real>::iterator beg,
			std::vector<real>::iterator end)
		{ v.erase(beg, end); }
		void insert(std::vector<real>::iterator beg,
			real val){v.insert(beg, val);}
		void insert(std::vector<real>::iterator beg, 
			std::vector<real>::iterator end,
			std::vector<real>::iterator last)
		{ v.insert(beg,end,last); }
		void push_back(real x) { v.push_back(x); };
		void resize(int i) { v.resize(i, 0); };
		void resize(int i, real j) { v.resize(i, j); };
		size_t size() { return v.size(); }

		// Overloaded operators.
		real operator [](int i) const { return v[i]; }
		real& operator [](int i) { return v[i]; }
		bool operator==(const Vector& rhs) { return v == rhs.v; }
		bool operator!=(const Vector& rhs) { return v != rhs.v; }
		Vector& operator*=(Vector& rhs);
		Vector& operator*=(real x);
		Vector& operator/=(Vector& rhs);
		Vector& operator/=(real x);
		Vector& operator+=(Vector& rhs);
		Vector& operator+=(real x);
		Vector& operator-=(Vector& rhs);
		Vector& operator-=(real x);
		friend Vector operator+(Vector& lhs, Vector& rhs);
		friend Vector operator+(Vector& lhs, real x);
		friend Vector operator+(real x, Vector& rhs);
		friend Vector operator-(Vector& lhs, Vector& rhs);
		friend Vector operator-(Vector& lhs, real x);
		friend Vector operator-(real x, Vector& rhs);
		friend Vector operator*(Vector& lhs, Vector& rhs);
		friend Vector operator*(Vector& lhs, real x);
		friend Vector operator*(real x, Vector& rhs);
		friend Vector operator/(Vector& lhs, Vector& rhs);
		friend Vector operator/(Vector& lhs, real x);
		friend Vector operator/(real x, Vector& rhs);

		// Unique class methods.
		std::vector<real>* get() { return &v; }
		real length();
		real mean();
		void removeNullValues();
		real populationVariance();
		real populationStandardDeviation();
		real sampleVariance();
		real sampleStandardDeviation();
		real RMS();
		real sum();
		real sumexp(real p);
		virtual std::string to_string(int precision = 4);

	private:
		// Wrapped std::vector containing real values.
		std::vector<real> v;

		// Helper functions to reduce reusing code for operator overloads.
		Vector add(Vector& a, Vector& b);
		Vector add(Vector& a, real  b);
		Vector add(real b, Vector& a);
		Vector subtract(Vector& a, Vector& b);
		Vector subtract(Vector& a, real  b);
		Vector subtract(real b, Vector& a);
		Vector multiply(Vector& a, Vector& b);
		Vector multiply(Vector& a, real b);
		Vector multiply(real b, Vector& a);
		Vector divide(Vector& p, Vector& q);
		Vector divide(Vector& p2, real q);
		Vector divide(real q, Vector& p2);
	};

	real dot(Vector v1, Vector v2);
	Vector ceil(Vector v);
	Vector clamp(Vector v1, real min_, real max_);
	Vector cross(Vector v1, Vector v2);
	Vector floor(Vector v);
	Vector fract(Vector v);
	real length(Vector v);
	Vector mix(Vector a, Vector b, real t);
	Vector normalize(Vector v);
	Vector round(Vector v);
}

namespace std
{
	MathParser::Vector pow(MathParser::Vector v, int n);
}