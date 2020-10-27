#pragma once
#include "Types.h"
#include <vector>

namespace MathParser
{
	class ComplexNumber;
	class Function;
	class Matrix;
	class Polynomial;
	class Vector;

	// Tables of weights for numerical integration schemes.
	extern const highp GaussLegendreTable[100];
	extern const highp GaussHermiteTable[150];
	extern const highp GaussLaguerreTable[256];
	extern const highp GaussLaguerreTableWithAlpha[300];
	extern const highp GaussKronodTable[198];
	extern const highp GaussKronodTableCorrespondingGaussTable[98];

	// Functions for integration/differentiation schemes.
	real nthDerivative1D(real D, real x, Function f, highp epsilon = 0.000001);
	real nthPartialDerivative(real D, int wrt, Vector x, Function f,highp epsilon = 0.000001);
	real integrateSimpsonsRule(real a, real b, Function f, real subdivisions = 2000.0);
	real integrateGaussLegendreQuadrature(real a, real b, Function f);
	real integrateGaussLegendreQuadrature(real a, real b, Function f, Vector vals, int pos);
	real integrateGaussHermiteQuadrature(Function f);
	real integrateGaussLaguerreQuadrature(Function f);
	real integrateGaussLaguerreQuadrature(Function f, Vector vals, int pos);
	real integrateGaussKronodQuadrature(real a, real b, Function f);
	real integrateGaussKronodQuadrature(real a, real b, Function f, Vector vals, int pos);
}