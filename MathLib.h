/**
*   Math library for use with the math parser.
*/

#pragma once
#include "Types.h"
#include <vector>

namespace MathParser
{
	int getMantissaLength(real x);
	highpUint numberOfDivisors(highpUint n);
	std::vector<highpUint> divisors(highpUint n);

	// Note: this method uses pre-calculated tables of primes for the first 3600 primes.
	// Beyond that, a more complex method is used for calculation.
	bool isPrime(highpUint n);
	int prime(int n);
	int sumPrimes(int n);
	int productPrimes(int n);
	int primeCountingFunction(int n);

	int AckermannsFunction(int m, int n);
	real BernoulliNumber(int n);
	real BernoulliNumber(int m, int n);
	highpInt CollatzConjecture(highpInt n);
	highpUint doubleFactorial(highpInt n);
	real factorial(real n);
	real FibonacciNumber(int n);
	bool isSquarefree(int n);
	real logarithm(real n, real base);
	int permutation(int n, int r);
	int combination(int n, int r);
	int MoebiusFunction(int n);
	real fallingFactorial(real n, real r);
	real risingFactorial(real n, real r);
	real digamma(real z);
	real gamma(real z);
	real lowerIncompleteGamma(real s, real x);
	real upperIncompleteGamma(real s, real x);
	real regularizedGamma(real s, real x);
	real logit(real x);
	real sigmoid(real n);
	real LambertW(real z);
	real WrightOmegaFunction(real x);
	highpUint BellNumber(real n);
	int StirlingNumber1stKind(int n, int k);
	int UnsignedStirlingNumber1stKind(int n, int k);
	real StirlingNumber2ndKind(real n, real k);
	int EulersTotientFunction(int n);
	real BesselFunction2ndKind(real order, real x);

	//form: (a;q)_n = Product k=[0,n-1] of (1-aq^k).
	real qPochhammerFactorial(real a, real q, int n);

	int sgn(real x);
	real step(real t);
	real rect(real x);
	real boxcar(real start, real stop, real amplitude, real t);
	real tri(real x);
	real sinc(real x);
	real squarewave(real t);
}