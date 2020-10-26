#pragma once
#include "Types.h"

namespace MathParser
{
	class ComplexNumber;
	class Function;
	class Matrix;
	class Polynomial;
	class Vector;

	Polynomial linearLeastSquaresRegression(Matrix A);
	Polynomial linearLeastSquaresRegression(Vector v1, Vector v2);
	real covariance(Vector v1, Vector v2);
	real normalCDF(real z);
	real normalCDFApproximation(real z);
	real normalPDF(real z);
	real CauchyPDF(real x, real x0, real gamma);
	real CauchyCDF(real x, real x0, real gamma);
	real studentsCDF(real t, real n, real epsilon=0.0000001);
	real studentsPDF(real x, real v);
	real chiSquareCDF(real k, real x);
	real ChiSquareTest(Matrix A);
	real lognormalDistributionPDF(real x, real mean, real variance);
	real lognormalDistributionCDF(real x, real mean, real variance);
	real FDistributionApproximation(real x, real v1, real v2);
	real GompertzPDF(real x, real n, real b);
	real GompertzCDF(real x, real n, real b);
	real correlationCoefficient(Vector a, Vector b);
	real coefficientOfDetermination(Vector x, Vector y);
	Matrix correlationMatrix(Matrix A);
	Matrix covarianceMatrix(Matrix M);
	Vector JarqueBeraTest(Matrix M, real criticalValue);
	Vector JarqueBeraSkewness(Matrix M);
	Vector JarqueBeraKurtosis(Matrix M);
	Vector DAgostinoTest(Matrix M, real criticalValue);
	Vector DAgostinoSkewness(Matrix M);
	Vector DAgostinoKurtosis(Matrix M);
	Vector AndersonDarlingTest(Matrix M);
	Vector AndersonDarlingTest(Vector X);
	real RMS(int T1, int T2, Function f);
}