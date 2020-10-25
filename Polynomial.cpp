#include "Polynomial.h"
#include "StringUtils.h"
#include "MathLib.h"
#include <limits>
#include <ctime>
#include <iomanip>
#include <cmath>

namespace MathParser
{
	Polynomial::Polynomial(std::vector<real> coef)
	{
		coefficient = coef;
	}

	int Polynomial::getLength(real x)
	{
		real n = abs(x);
		int length = 0;
		if (n < 1)
			return 0;
		while (n > 1)
		{
			n /= 10;
			++length;
		}
		return length;
	}

	void Polynomial::randomize() 
	{
		srand(clock());
		int terms = coefficient.size();
		coefficient.clear();
		for (int i = 0; i < terms; i++) {
			real n = ((real )rand() / 10000) * pow(-1, rand() % 2);
			coefficient.push_back(n);
		}
	}

	real Polynomial::largestCoefficient() 
	{
		real biggest = 0;
		for (int i = 0; i < coefficient.size(); ++i) {
			if (fabsf(coefficient[i]) > biggest) { biggest = abs(coefficient[i]); }
		}
		return biggest;
	}
	
	real Polynomial::evaluate(real n, int precision) 
	{// Uses Horner's method for optimal speed.
		real result = 0.0;
		for (int idx = coefficient.size() - 1; idx >= 0; idx--) 		
			result = fma(result, n, coefficient[idx]);		
		return result;
	}

	complex Polynomial::evaluate(complex z) 
	{
		std::vector<real> v = coefficient;
		complex s(0, 0);
		for (std::vector<real>::const_reverse_iterator i = v.rbegin(); i != v.rend(); i++)
			s = (s * z) + *i;		
		return s;
	}

	
	real Polynomial::NewtonsMethod(real x) 
	{// Given some root guess, this method will refine it until it is much more accurate.
		real x1;
		real f, dfdx;
		const real epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = derivative();

		x1 = x;//initial guess

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = evaluate(x);
			dfdx = deriv.evaluate(x);
			x1 = x - f / dfdx;
		}
		
		if (abs(x - x1) < epsilon) 
			return x1;

		return 0;
	}

	bool Polynomial::isIntegerValued() {
		for (int i = 0; i < coefficient.size(); ++i) {
			if (coefficient[i] != floor(coefficient[i])) { return false; }
		}
		return true;
	}

	bool Polynomial::isMonic() {
		if (coefficient.size() > 0 && coefficient[coefficient.size() - 1] != 1) { return false; }
		return true;
	}

	Polynomial Polynomial::makeMonic() {
		if (coefficient[coefficient.size() - 1] != 1) {
			real a = 1 / coefficient[coefficient.size() - 1];
			for (int i = 0; i < coefficient.size(); ++i) { coefficient[i] *= a; }
		}
		return *this;
	}

	Polynomial Polynomial::factorOutBinomial(Polynomial p) 
	{
		real n = p.coefficient[0] - .1; //get -n from the input binomial (x + n)
		if (p.coefficient.size() > 2) { return Polynomial(); }//can't factor out larger polynomials

		//synthetic division w/ adjustment for imprecise 
		//================================================
		std::vector<real> vec;
		int counter = 0;
		real remainder = 99999;
		real delta = 0.001;
		int sign = 1;
		while (std::abs(remainder) > 0.001 && counter < 500) 
		{
			vec.resize(coefficient.size() - 1);

			// Drop down first coefficient of synth div.
			vec[coefficient.size() - 2] = coefficient[coefficient.size() - 1];
			
			// Do division.
			for (int i = coefficient.size() - 2; i >= 0; --i) {
				if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
				else { vec[i - 1] = coefficient[i] + (vec[i] * n); }
			}
			
			++counter;
			
			if (counter <= 1) 
			{ 
				sign = sgn(remainder); 
			}

			if (counter > 1) 
			{ 
				if (sign != sgn(remainder)) 
				{ 
					n -= delta; delta /= 10; 
				} 
			}

			if (std::abs(remainder) < 0.00001)			
				return Polynomial(vec);			
			if (std::abs(remainder) >= 0.00001) 
				n += delta;
		}
		return Polynomial(vec);
	}

	Polynomial Polynomial::factorOutBinomial(real n) {//factors out a binomial by synthetic division
		std::vector<real > cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return factorOutBinomial(cof);
	}

	real Polynomial::syntheticDivision(Polynomial p) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the remainder
		real n = p.coefficient[0]; //get -n from the input binomial (x + n)
		
		// For now, can't factor out larger polynomials.
		if (p.coefficient.size() > 2)
			return 0;

		//synthetic division w/ adjustment for imprecise n
		//================================================
		int counter = 0;
		real remainder;
		real delta = 0.001;
		int sign = 1;
		std::vector<real> vec;
		vec.resize(coefficient.size() - 1);

		// Drop down first coefficient of synth div.
		vec[coefficient.size() - 2] = coefficient[coefficient.size() - 1];

		// Do division, iterate from largest to smallest exponent in polynomial.
		for (int i = coefficient.size() - 2; i >= 0; --i) {
			if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
			else { vec[i - 1] = coefficient[i] + (vec[i] * n); }
		}

		return remainder;
	}

	real Polynomial::syntheticDivision(real n) {//factors out a binomial by synthetic division, return remainder
		std::vector<real > cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return syntheticDivision(cof);
	}

	complex Polynomial::complexSyntheticDivision(complex z) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the remainder
		int counter = 0;
		complex remainder;
		std::vector<complex> vec;
		vec.resize(coefficient.size() - 1);

		vec[coefficient.size() - 2] = complex(coefficient[coefficient.size() - 1], 0);//drop down first coefficient of synth div
		for (int i = coefficient.size() - 2; i >= 0; --i) {
			if (i == 0) 
				remainder = coefficient[i] + (z * vec[i]);
			else
				vec[i - 1] = coefficient[i] + (vec[i] * z);			
		}

		return remainder;
	}

	bool Polynomial::checkLowBound(Polynomial p) {//factors out a binomial by synthetic division where binomial p = (x + ...) and returns the true is signs alternate
		real n = p.coefficient[0]; //get n from the input binomial (x + n)
		if (p.coefficient.size() > 2) { return 0; }//can't factor out larger polynomials

									   //synthetic division
									   //==================
		real * vec = NULL;
		int counter = 0;
		real remainder;
		real delta = 0.001;
		int signflips = 0;
		int sign = 0;

		vec = new real [coefficient.size() - 1];
		vec[coefficient.size() - 2] = coefficient[coefficient.size() - 1];//drop down first coefficient of synth div
		sign = sgn(vec[coefficient.size() - 2]);
		for (int i = coefficient.size() - 2; i >= 0; --i) {
			if (i == 0) { remainder = coefficient[i] + (vec[i] * n); }
			else {
				vec[i - 1] = coefficient[i] + (vec[i] * n);
				int sn = sgn(vec[i - 1]);
				if (sign != sn) {
					++signflips;
					sign = sn;
				}
			}
		}
		if (sgn(remainder) != sign) { ++sign; }
		if (signflips >= coefficient.size() - 2) { return true; }
		return false;
	}

	bool Polynomial::checkLowBound(real n) {//factors out a binomial by synthetic division, return remainder
		std::vector<real> cof;
		cof.push_back(n);
		cof.push_back(1);
		Polynomial p(cof);
		return checkLowBound(cof);
	}

	Polynomial Polynomial::derivative() {
		if (!coefficient.size() <= 1) { return Polynomial(); }
		std::vector<real> coef1;
		coef1.resize(coefficient.size()-1);
		for (int i = 0; i < coefficient.size() - 1; i++) {
			coef1[i] = coefficient[i + 1] * (i + 1);
		}
		return Polynomial(coef1);
	}

	Polynomial Polynomial::integral() {
		if (!coefficient.size() <= 1) { return Polynomial(); }
		std::vector<real> coef1;
		coef1.resize(coefficient.size() + 1);
		for (int i = 0; i < coefficient.size()+1; i++) {
			coef1[i] = coefficient[i + 1] * (i + 1);
		}
		return Polynomial(coef1);
	}

	std::vector<real> Polynomial::getRootBounds() {
		Polynomial temp(coefficient);
		if (temp.isMonic() == false) { temp.makeMonic(); }
		real M = temp.largestCoefficient() + 1;//by Cauchy's theorem
		std::vector<real > answer;
		int lb = 0;
		int hb = 0;

		// Now, test for smaller bounds by using synthetic division:
		//check lower bound
		for (real i = -M; i < M; ++i) {
			if (checkLowBound(i)) 
				lb = i;
			if (!checkLowBound(i)) //once you've passed beyond the bound, break
				i = M; 
		}
		if (lb != 0) { answer.push_back(lb); }
		else { answer.push_back(-M); }

		//check upper bound
		for (real i = 0; i < M; i++) {
			if (std::abs(syntheticDivision(i)) < 0.01) { hb = i; i = M; }
		}
		if (hb != 0) { answer.push_back(hb); }
		else { answer.push_back(M); }

		return answer;
	}

	real Polynomial::root(Polynomial p, real lower_bound, real upper_bound) {//bounded search for a root
		//given a set of bounds for a polynomial, this method will find the first root
		if (p.coefficient.size() <= 1) { return 0; }
		if (p.coefficient.size() == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.coefficient.size() == 3) {
			real part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		real TOL = 0.0001;    // tolerance
		real MAX_ITER = 1000; // maximum number of iterations
		real a = lower_bound;
		real b = upper_bound;
		real fa = p.evaluate(a);   // calculated now to save function calls
		real fb = p.evaluate(b);   // calculated now to save function calls
		real fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		real c = a;           // c now equals the largest magnitude of the lower and upper bounds
		real fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		real s = 0;           // Our Root that will be returned
		real d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

		  /*	Condition statement:
		  -------------------------------------------------------
		  (condition 1) s is not between  (3a+b)/4  and b or
		  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
		  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
		  (condition 4) (mflag is set and |b−c| < |TOL|) or
		  (condition 5) (mflag is false and |c−d| < |TOL|) */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b) * 0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.evaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		real x, x1; //user input value
		real f, dfdx; //function and its derivative 
		const real epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.evaluate(x); //Defines the function 
			dfdx = deriv.evaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	std::vector<real> Polynomial::realRoots() {
		Polynomial p(coefficient);//Make a copy of this polynomial for alteration.
		
		std::vector<real> candidates;
		if (p.coefficient.size() <= 1) { return candidates; }
		if (p.coefficient.size() == 2) {
			candidates.push_back(-p.coefficient[0] / p.coefficient[1]);
			return candidates;
		}
		if (p.coefficient.size() == 3) {
			real part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 >= 0) {
				real part1 = (-p.coefficient[1]) / (2 * p.coefficient[2]);
				part2 = sqrt(part2) / (2 * p.coefficient[2]);
				candidates.push_back(part1 + part2);
				candidates.push_back(part1 - part2);
			}
			return candidates;
		}
		//if the polynomial is longer then 2nd order...
		std::vector<real > bounds = p.getRootBounds();
		real lowbound = bounds[0];
		real highbound = bounds[1];

		//Since 0 is used as a null value for roots in the rootfinding algorithms,
		//we must determine whether or not it is a root outside of these algorithms.
		if (p.evaluate(0) == 0) { candidates.push_back(0); }

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		Polynomial temp = p;
		real lastroot = 0;
		real delta = 0.01;
		if (getLength(p.largestCoefficient()) > 3) { delta *= 10 * getLength(p.largestCoefficient()); }
		while (lowbound < highbound) {
			real x = root(temp, lowbound, lowbound + delta);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						bounds = temp.getRootBounds();
						real lowbound = bounds[0];
						real highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	std::vector<real> Polynomial::getComplexRootBounds() {
		Polynomial temp(coefficient);
		if (temp.isMonic() == false) { temp.makeMonic(); }
		real M = temp.largestCoefficient() + 1;//by Cauchy's theorem
		std::vector<real > answer;
		real Hbound = 0.5 * ((std::abs(coefficient[coefficient.size() - 2]) - 1) + 
			sqrt(pow(coefficient[coefficient.size() - 2] - 1, 2) + (4 * M)));//Sun and Hsieh bound
		answer.push_back(-Hbound);
		answer.push_back(Hbound);

		std::vector<real > co;
		co.push_back(1);
		co.push_back(2 - std::abs(coefficient[coefficient.size() - 2]));
		co.push_back(1 - std::abs(coefficient[coefficient.size() - 2]) - 
			std::abs(coefficient[coefficient.size() - 3]));
		co.push_back(-M);
		Polynomial P(co);
		co = P.realRoots();
		if (co.size() > 1) {
			answer[0] = co[0];
			answer[1] = co[1];
		}
		return answer;
	}

	std::string Polynomial::to_string(int precision) {
		std::stringstream s;

		int lastnonzero = 0;
		for (int i = coefficient.size() - 1; i >= 0; i--) { 
			if (coefficient[i] != 0) { lastnonzero = i; } }

		for (int i = coefficient.size() - 1; i >= 0; i--) {
			if (coefficient[i] != 0) {
				if (coefficient[i] == floor(coefficient[i])) { precision = 0; }
				if (coefficient[i] != floor(coefficient[i])) { precision = precision; }

				if (i == 0) {
					if (coefficient[i] >= 0) {
						s << " + " << std::setprecision(precision) << coefficient[i];
					}
					if (coefficient[i] < 0) {
						std::string str = " - ";
						if (i == coefficient.size() - 1) { str = "-"; }
						s << str << std::setprecision(precision) << std::abs(coefficient[i]);
					}
				}
				if (i == 1) {
					if (std::abs(coefficient[i]) != 1) {
						if (coefficient[i] >= 0) {
							s << " + " << std::setprecision(precision) << coefficient[i] << "x";
						}
						if (coefficient[i] < 0) {
							std::string str = " - ";
							if (i == coefficient.size() - 1) { str = "-"; }
							s << str << std::setprecision(precision) << std::abs(coefficient[i]) << "x";
						}
					}
					if (coefficient[i] == 1) {
						s << " + " << "x";
					}
					if (coefficient[i] == -1) {
						std::string str = " - ";
						if (i == coefficient.size() - 1) { str = " - "; }
						s << str << "x";
					}
				}
				if (i > 1) {
					if (coefficient[i] != 1 && coefficient[i] != -1) {

						if (coefficient[i] >= 0) {
							s << " + " << std::setprecision(precision) << coefficient[i] << "x"
								<< "^" << i;
						}
						if (coefficient[i] < 0) {
							std::string str = " - ";
							if (i == coefficient.size() - 1) { str = "-"; }
							s << str << std::setprecision(precision) << std::abs(coefficient[i]) 
								<< "x^" << i;
						}

					}
					if (coefficient[i] == 1) {
						s << " + " << "x^" << i;
					}
					if (coefficient[i] == -1) {
						std::string str = " - ";
						if (i == coefficient.size() - 1) { str = "-"; }
						s << str << "x^" << i;
					}
				}
			}
		}
		std::string temp = s.str();
		while (!std::isdigit(temp[0]) && temp[0] != 'x' && temp[0] != '-')
			temp = temp.substr(1);
		return temp;
	}

	real Polynomial::findRoot() {
		//Make copy of this polynomial for alteration.
		Polynomial p(coefficient);
		
		//given a set of bounds for a polynomial, this method will find the first root
		if (p.coefficient.size() <= 1) { return 0; }
		if (p.coefficient.size() == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.coefficient.size() == 3) {
			real part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		real bound = 1 + (fabsf(p.largestCoefficient()));
		real TOL = 0.0001;    // tolerance
		real MAX_ITER = 1000; // maximum number of iterations
		real a = -bound; //lower_bound;
		real b = bound; //upper_bound;
		real fa = p.evaluate(a);   // calculated now to save function calls
		real fb = p.evaluate(b);   // calculated now to save function calls
		real fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		real c = a;           // c now equals the largest magnitude of the lower and upper bounds
		real fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		real s = 0;           // Our Root that will be returned
		real d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

													  /*	  (condition 1) s is not between  (3a+b)/4  and b or
													  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
													  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
													  (condition 4) (mflag is set and |b−c| < |TOL|) or
													  (condition 5) (mflag is false and |c−d| < |TOL|)
													  */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b) * 0.5;
				mflag = true;
			}
			else { mflag = false; }

			fs = p.evaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		real x, x1; //user input value
		real f, dfdx; //function and its derivative 
		const real epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.evaluate(x); //Defines the function 
			dfdx = deriv.evaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	real Polynomial::findRoot(real lower_bound, real upper_bound) {//bounded search for a root
																		   //given a set of bounds for a polynomial, this method will find the first root
		Polynomial p(coefficient);

		if (p.coefficient.size() <= 1) { return 0; }
		if (p.coefficient.size() == 2) { return (-p.coefficient[0] / p.coefficient[1]); }
		if (p.coefficient.size() == 3) {
			real part1 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part1 >= 0) { return ((-p.coefficient[1] + sqrt(part1)) / (2 * p.coefficient[2])); }
			else { return 0; }
		}

		//else, p is of order 3 or higher....
		real TOL = 0.0001;    // tolerance
		real MAX_ITER = 1000; // maximum number of iterations
		real a = lower_bound;
		real b = upper_bound;
		real fa = p.evaluate(a);   // calculated now to save function calls
		real fb = p.evaluate(b);   // calculated now to save function calls
		real fs = 0;      // initialize 

		if (!(fa * fb < 0)) { return 0; }

		if (fabsf(fa) < fabsf(b)) { // if magnitude of f(lower_bound) is less than magnitude of f(upper_bound)
			std::swap(a, b);
			std::swap(fa, fb);
		}

		real c = a;           // c now equals the largest magnitude of the lower and upper bounds
		real fc = fa;         // precompute function evalutation for point c by assigning it the same value as fa
		bool mflag = true;      // boolean flag used to evaluate if statement later on
		real s = 0;           // Our Root that will be returned
		real d = 0;           // Only used if mflag is unset (mflag == false)

		for (unsigned int iter = 1; iter < MAX_ITER; ++iter)
		{
			// stop if converged on root or error is less than tolerance
			if (std::abs(b - a) < TOL) {
				goto NarrowResult;
			} // end if

			if (fa != fc && fb != fc) {
				// use inverse quadratic interopolation
				s = (a * fb * fc / ((fa - fb) * (fa - fc)))
					+ (b * fa * fc / ((fb - fa) * (fb - fc)))
					+ (c * fa * fb / ((fc - fa) * (fc - fb)));
			}
			else { s = b - fb * (b - a) / (fb - fa); }//secant method

		  /* Condition statement:
		   ----------------------
		  (condition 1) s is not between  (3a+b)/4  and b or
		  (condition 2) (mflag is true and |s−b| ≥ |b−c|/2) or
		  (condition 3) (mflag is false and |s−b| ≥ |c−d|/2) or
		  (condition 4) (mflag is set and |b−c| < |TOL|) or
		  (condition 5) (mflag is false and |c−d| < |TOL|)  */
			if (((s < (3 * a + b) * 0.25) || (s > b)) ||
				(mflag && (std::abs(s - b) >= (std::abs(b - c) * 0.5))) ||
				(!mflag && (std::abs(s - b) >= (std::abs(c - d) * 0.5))) ||
				(mflag && (std::abs(b - c) < TOL)) ||
				(!mflag && (std::abs(c - d) < TOL)))
			{
				// bisection method
				s = (a + b) * 0.5;
				mflag = true;
			}
			else
				mflag = false;

			fs = p.evaluate(s);  // calculate fs
			d = c;      // first time d is being used (wasnt used on first iteration because mflag was set)
			c = b;      // set c equal to upper bound
			fc = fb;    // set f(c) = f(b)

			if (fa * fs < 0) {// fa and fs have opposite signs
				b = s;
				fb = fs;    // set f(b) = f(s)
			}
			else {
				a = s;
				fa = fs;    // set f(a) = f(s)
			}

			if (std::abs(fa) < std::abs(fb)) {// if magnitude of fa is less than magnitude of fb
				std::swap(a, b);     // swap a and b
				std::swap(fa, fb);   // make sure f(a) and f(b) are correct after swap
			}
		}
		return 0;

	NarrowResult:
		//use root as guess for Newton's method now
		real x, x1; //user input value
		real f, dfdx; //function and its derivative 
		const real epsilon = 0.0000001; //tolerance
		const int n = 15; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = s;//first guess is s

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.evaluate(x); //Defines the function 
			dfdx = deriv.evaluate(x); //Derivative of the function
			x1 = x - f / dfdx; //Newton Method formula
		}
		if (abs(x - x1) < epsilon) { //checks if guesses are within a certain tolerance value 
			return x1;
		}

		return 0;
	}

	std::vector<real> Polynomial::findRealRoots() {
		Polynomial p(coefficient);

		std::vector<real> candidates;
		if (p.coefficient.size() <= 1) { return candidates; }
		if (p.coefficient.size() == 2) {
			candidates.push_back(-p.coefficient[0] / p.coefficient[1]);
			return candidates;
		}
		if (p.coefficient.size() == 3) {
			real part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 >= 0) {
				real part1 = (-p.coefficient[1]) / (2 * p.coefficient[2]);
				part2 = sqrt(part2) / (2 * p.coefficient[2]);
				candidates.push_back(part1 + part2);
				candidates.push_back(part1 - part2);
			}
			return candidates;
		}
		//if the polynomial is longer then 2nd order...
		std::vector<real > bounds = p.getRootBounds();
		real lowbound = bounds[0];
		real highbound = bounds[1];

		// 0 is used as a null value for roots in the rootfinding algorithms.
		if (p.evaluate(0) == 0) { candidates.push_back(0); }

		//now sweep through the bounds of the polynomial to look for roots using the
		//mixed Brent/Newton search method in findRoot()
		Polynomial temp = p;
		real lastroot = 0;
		real delta = 0.01;
		if (getLength(p.largestCoefficient()) > 3) { delta *= 10 * getLength(p.largestCoefficient()); }
		while (lowbound < highbound) {
			real x = temp.findRoot(lowbound, lowbound + delta);
			if (x == 0) { lowbound += delta; }//using 0 as a null return value here
			else {
				if (x != lastroot) { //prevent same root from continually showing up
					candidates.push_back(x);
					lastroot = x;
					lowbound = x + delta;
					if (x == floor(x)) {
						temp = temp.factorOutBinomial(x); //perform a shift to reduce the polynomial
						bounds = temp.getRootBounds();
						real lowbound = bounds[0];
						real highbound = bounds[1];
					}//if x is an integer, you can factor it out of the polynomial
				}
				else { lowbound += delta; }
			}
		}
		return candidates;
	}

	complex Polynomial::NewtonsMethodComplex(Polynomial p, complex z) {
		//use root as guess for Newton's method now
		complex x, x1; //user input value
		complex f, dfdx; //function and its derivative 
		const real epsilon = 0.00000001; //tolerance
		const int n = 1000; //# steps
		Polynomial deriv = p.derivative();//derivative of polynomial

		x1 = z;//initial guess.

		for (int i = 2; i <= n; ++i) {
			x = x1;
			f = p.evaluate(x); //Defines the function 
			dfdx = deriv.evaluate(x); //Derivative of the function
			x1 = x - (f / dfdx); //Newton Method formula = x1 = x - f/dfdx
		}

		if (x.imag() != 0 && (x - x1).real() < epsilon && (x - x1).imag() < epsilon) { //checks if guesses are within a certain tolerance value AND complex (not real)
			return x1;
		}
		return complex();//if not close enough, return 0 value
	}
	

	std::vector<complex> Polynomial::roots() {
		Polynomial p(coefficient);

		std::vector<real> realroots = p.findRealRoots();
		std::vector<complex> complexRoots;
		for (int i = 0; i < realroots.size(); ++i)
			complexRoots.push_back(complex(realroots[i], 0));		

		std::vector<real > bounds = p.getComplexRootBounds();
		real lowbound = bounds[0];
		real highbound = bounds[1];
		real lowboundC = bounds[0];

		for (int i = 0; i < realroots.size(); ++i) {
			if (std::abs(int(p.coefficient.size()) - int(realroots.size())) <= 3) {
				// Use the shift method only if the real roots reduce it down to an 
				// easily calculable size (order 2 or less),
				// otherwise the results will be inaccurate
				if (p.coefficient.size() <= 3) { i = realroots.size(); break; }
				if (p.coefficient.size() > 3) { p = p.factorOutBinomial(realroots[i]); }
			}
		}
		if (p.coefficient.size() < 3) { return complexRoots; }

		if (p.coefficient.size() == 3) {
			real part2 = pow(p.coefficient[1], 2) - (4 * p.coefficient[0] * p.coefficient[2]);
			if (part2 < 0) {
				part2 = part2 / (2 * p.coefficient[2]);
				real part1 = (-p.coefficient[1] / (2 * p.coefficient[2]));
				complexRoots.push_back(complex(part1, part2));
				complexRoots.push_back(complex(part1, -part2));
			}
			return complexRoots;
		}

		if (p.coefficient.size() > 3) {

			//start finding roots with Durand-Kerner method
			int iterations = 10000;
			complex z = complex(lowbound, lowboundC);
			int size = sizeof(z);
			std::vector<complex> R;
			for (int i = 0; i < p.coefficient.size(); i++) { R.push_back(std::pow(z,i)); }

			for (int i = 0; i < iterations; i++) {
				for (int j = 0; j < p.coefficient.size(); j++) {
					complex B = p.evaluate(R[j]);
					for (int k = 0; k < p.coefficient.size(); k++) {
						if (k != j) { B /= (R[j] - R[k]); }
					}
					R[j] -= B;
				}
			}

			for (int i = 0; i < R.size(); ++i) { //now filter out all the bad results
				if (R[i].imag() != 0) { //only complex roots accepted
					R[i] = NewtonsMethodComplex(p, R[i]);
					complex temp = p.evaluate(R[i]);
					real a = std::abs(temp.real());
					real b = std::abs(temp.imag());
					if (a < 0.1 && b < 0.1) {
						a = std::abs(R[i].real());
						b = std::abs(R[i].imag());
						bool isSimilar = false;
						for (int i = 0; i < complexRoots.size(); ++i) {
							real x = std::abs(complexRoots[i].real());
							real y = std::abs(complexRoots[i].imag());
							if (std::abs(a - x) < 0.001 && std::abs(b - y) < 0.001) { isSimilar = true; }
						}
						if (isSimilar == false) { //if this is indeed a new root, save the root and its conjugate
							complexRoots.push_back(R[i]);
							complexRoots.push_back(complex(R[i].real(), R[i].imag() * -1));
						}
					}
				}
			}
			R.clear();
		}
		return complexRoots;
	}
}