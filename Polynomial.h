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
	class Polynomial : public BaseObject
	{
	public:
		Polynomial() = default;
		Polynomial(std::vector<real> coef);

		virtual real evaluate(real n, int precision = 4);
		complex evaluate(complex z);
		std::vector<real> realRoots();
		virtual std::string to_string(int precision=4);
		std::vector<complex> roots();

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

	protected:
		// Array of per-exponent coefficient values for the polynomial.
		std::vector<real> coefficient;

	private:
		// For ease of testing, randomize this polynomial's terms and coefficients.
		void randomize();

		// Root-finding helper methods //
		// Factors out a binomial by synthetic division where binomial p = (x + ...).
		Polynomial factorOutBinomial(Polynomial p);

		// Find largest coefficient of polynomial.
		real largestCoefficient();

		// Helper function to get the length of the largest coefficient > 1.
		int getLength(real x);

		bool isIntegerValued();
		bool isMonic();
		Polynomial makeMonic();
		Polynomial derivative();
		Polynomial integral();
		real NewtonsMethod(real x);
		complex NewtonsMethodComplex(Polynomial p, complex z);
		Polynomial factorOutBinomial(real n);
		real syntheticDivision(Polynomial p);
		real syntheticDivision(real n);
		complex complexSyntheticDivision(complex z);
		bool checkLowBound(real n);
		bool checkLowBound(Polynomial p);
		std::vector<real> getRootBounds();
		real root(Polynomial p, real lower_bound, real upper_bound);
		std::vector<real> getComplexRootBounds();
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
		//Polynomial divide(real q, real p2);
	};	
}

namespace std
{
	// overload for standard 'pow()' function to use Polynomial class.
	MathParser::Polynomial pow(MathParser::Polynomial p, int b);
}