/**
*  \brief Polynomial class used for polynomials with positive, integer-valued
*		  exponents, represented as indices of the coefficient vector 
*		  (ie '2x^2 + x + 1' is represented as 'coefficient' array {2,1,1}')
*/

#pragma once
#include <vector>
#include "Types.h"
#include "BaseObject.h"

namespace MathParser
{
	class ComplexNumber;
	class Matrix;
	class Vector;

	class Polynomial : public BaseObject
	{
	public:
		Polynomial() {}
		Polynomial(std::vector<real> coef);
		Polynomial(Vector v);

		// Overloads of class functions for std::vector for brevity.
		void clear() { coefficient.clear(); }
		void push_back(real x) { coefficient.push_back(x); }
		size_t size() { return coefficient.size(); }

		// Overloaded operators.
		real operator [](int i) const { return coefficient[i]; }
		real& operator [](int i) { return coefficient[i]; }
		bool operator==(const Polynomial& rhs) { return coefficient == rhs.coefficient; }
		Polynomial& operator*=(Polynomial& rhs);
		Polynomial& operator*=(real x);
		Polynomial& operator/=(Polynomial& rhs);
		Polynomial& operator/=(real x);
		Polynomial& operator+=(Polynomial& rhs);
		Polynomial& operator+=(real x);
		Polynomial& operator-=(Polynomial& rhs);
		Polynomial& operator-=(real x);
		friend Polynomial operator+(Polynomial& lhs, Polynomial rhs);
		friend Polynomial operator+(Polynomial& lhs, real x);
		friend Polynomial operator+(real x, Polynomial& rhs);
		friend Polynomial operator-(Polynomial& lhs, Polynomial rhs);
		friend Polynomial operator-(Polynomial& lhs, real x);
		friend Polynomial operator-(real x, Polynomial& rhs);
		friend Polynomial operator*(Polynomial& lhs, Polynomial rhs);
		friend Polynomial operator*(Polynomial& lhs, real x);
		friend Polynomial operator*(real x, Polynomial& rhs);
		friend Polynomial operator/(Polynomial& lhs, Polynomial rhs);
		friend Polynomial operator/(Polynomial& lhs, real x);
		//friend Polynomial operator/(real x, Polynomial& rhs);

		Polynomial derivative();
		virtual real evaluate(real n, int precision = 4);
		ComplexNumber evaluate(ComplexNumber z);		
		Polynomial factorOutBinomial(Polynomial p);// Factors out a binomial by synthetic division.
		Polynomial factorOutBinomial(real n);
		std::vector<real> getComplexRootBounds();
		Polynomial integral();
		bool isIntegerValued();
		bool isMonic();
		real largestCoefficient();
		Polynomial makeMonic();
		ComplexNumber NewtonsMethodComplex(Polynomial p, ComplexNumber z);
		std::vector<real> realRoots();
		std::vector<ComplexNumber> roots();
		virtual std::string to_string(int precision = 4);


	protected:
		// Wrapped array of per-exponent coefficient values for the polynomial.
		std::vector<real> coefficient;

	private:
		// For ease of testing, randomize this polynomial's terms and coefficients.
		void randomize();

		// Root-finding helper methods //
		real NewtonsMethod(real x);
		real syntheticDivision(Polynomial p);
		real syntheticDivision(real n);
		ComplexNumber complexSyntheticDivision(ComplexNumber z);
		bool checkLowBound(real n);
		bool checkLowBound(Polynomial p);
		std::vector<real> getRootBounds();
		real root(Polynomial p, real lower_bound, real upper_bound);
		real findRoot();
		real findRoot(real lower_bound, real upper_bound);
		std::vector<real> findRealRoots();

		// Helper functions to reduce reusing code for operator overloads.
		Polynomial add(Polynomial a, Polynomial b);
		Polynomial add(Polynomial a, real  b);
		Polynomial add(real b, Polynomial a);
		Polynomial subtract(Polynomial a, Polynomial b);
		Polynomial subtract(Polynomial a, real  b);
		Polynomial subtract(real b, Polynomial a);
		Polynomial multiply(Polynomial a, Polynomial b);
		Polynomial multiply(Polynomial a, real b);
		Polynomial multiply(real b, Polynomial a);
		Polynomial divide(Polynomial p, Polynomial q);
		Polynomial divide(Polynomial p2, real q);
		//Polynomial divide(real q, Polynomial p2);
	};	
}

namespace std
{
	// overload for standard 'pow()' function to use Polynomial class.
	MathParser::Polynomial pow(MathParser::Polynomial p, int b);
}