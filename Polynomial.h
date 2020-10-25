/**
*  \brief Polynomial class used for polynomials with positive, integer-valued exponents, represented
*			as indices of the coefficient vector (ie '2x^2' is represented as 'coefficient[3]=2')
*
*/

#pragma once
#include <unordered_map>
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
	};

}