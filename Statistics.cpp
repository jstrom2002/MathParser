#include "Statistics.h"
#include "Matrix.h"
#include "Polynomial.h"
#include "Vector.h"
#include "MathLib.h"
#include "ComplexNumber.h"
#include "Function.h"
#include "Calculus.h"
#include <algorithm>

namespace MathParser
{
	Polynomial linearLeastSquaresRegression(Matrix A) 
	{   //returns linear polynomial function
		if (A.columns == 2 && A.rows > 0) 
			A = A.transpose();

		//only accepts an input 2 x n matrix: one row for x values, the other for y values
		if (A.rows != 2) 
			return Polynomial(); 
		
		real sum_x = A.row(0).sum();
		real sum_y = A.row(1).sum();
		real x_mean = sum_x / A.columns;
		real y_mean = sum_y / A.columns;
		real sum_x_sqr = pow(A.row(0).length(), 2);
		real sum_xy = dot(A.row(0), A.row(1));

		real m = (sum_xy - ((sum_x * sum_y) / A.columns)) / 
			(sum_x_sqr - (sum_x * sum_x) / A.columns);
		real b = y_mean - (x_mean * m);

		std::vector<real> coef;
		coef.push_back(b);
		coef.push_back(m);
		std::vector<real> expo;
		expo.push_back(0);
		expo.push_back(1);

		Polynomial Regression(coef);
		return Regression;
	}

	Polynomial linearLeastSquaresRegression(Vector v1, Vector v2) 
	{// Handle data vector size mismatch.
		if (v1.size() != v2.size()) 
			return Polynomial(); 

		Vector vnew = v1;
		for (int i = 0; i < v2.size(); ++i)
			vnew.push_back(v2[i]);		
		Vector v3(vnew);
		Matrix M(v3);
		return linearLeastSquaresRegression(M);
	}

	real covariance(Vector v1, Vector v2) 
	{//sample covariance
		int sz = v1.size();
		real mu1 = v1.mean();
		real mu2 = v2.mean();
		real answer = 0;
		for (int i = 0; i < sz; ++i)
			answer += (v1[i] - mu1) * (v2[i] - mu2);		
		return answer / (sz - 1);
	}

	real normalCDFApproximation(real z) 
	{
		if (z == 0) 
			return 0.5;
		else //Vasquez-Leal normal cdf approximation
			return (pow((exp((-358.0 * z / 23.0) + 
				(111.0 * atan((37.0 * z) / 294.0))) + 1.0), -1.0));
	}

	real normalCDF(real x)
	{
		if (x == 0) 
			return 0.5;
		bool isNeg = x < 0;
		x = std::abs(x);
		real result = 0;
		real partitionwidth = .000001;
		real f = 0;
		for (real b = 0; b < x; b += partitionwidth) 
		{
			f = normalPDF(b);
			result += f;
		}
		result *= partitionwidth;
		if (isNeg) 
			return 1 - result;
		return result + 0.5;
	}

	real normalPDF(real z)
	{ 
		return ((1.0 / (2.0 * PI)) * exp(-0.5 * (z * z))); 
	}

	real studentsPDF(real x, real v) {
		if (v > 30) 
			return normalPDF(x);
		real temp = (v + 1) * 0.5;
		return ((tgamma(temp) / (tgamma(v * 0.5) * sqrt(v * PI))) * 
			(pow(v, temp) / pow((x * x) + v, temp)));
	}

	real studentsCDF(real t, real v, real epsilon) 
	{
		if (!t)
			return 0.5;
		real x = (v / ((t * t) + v));
		return (1 - (0.5 * incompleteBetaFunction(x, 0.5 * v, 0.5, epsilon)));
	}

	real Gaussian1D(real x, real o) {
		return exp(-0.5 * pow(x / o, 2)) / (o * sqrt(2 * PI));
	}

	real Gaussian1D(real x, real u, real o) {
		return exp(-0.5 * pow((x - u) / o, 2)) / (o * sqrt(2 * PI));
	}

	real Gaussian2D(real x, real y, real o) {
		return (exp(-1.0 * ((pow(x, 2) + pow(y, 2)) / (2 * pow(o, 2)))) / ((2 * PI) * pow(o, 2)));
	}

	real Gaussian2D(real x, real y, real u1, real u2, real o) {
		return (exp(-1.0 * ((pow(x - u1, 2) + pow(y - u2, 2)) / (2 * pow(o, 2)))) / ((2 * PI) * pow(o, 2)));
	}

	real GaussianND(std::vector<real> x, std::vector<real> o) {
		if (x.size() != o.size()) { return 0; }
		real val = 1;
		for (int i = 0; i < x.size(); ++i) {
			val *= Gaussian1D(x[i], o[i]);
		}
		return val;
	}

	real GompertzPDF(real x, real n, real b)
	{
		if (x < 0 || n <= 0 || b <= 0)
			return 0;
		return b * n * exp(n) * exp(b * x) * exp(-n * exp(b * x));
	}

	real GompertzCDF(real x, real n, real b)
	{
		if (x < 0 || n <= 0 || b <= 0)
			return 0;
		return 1 - exp(-n * (exp(b * x) - 1));
	}

	real CauchyPDF(real x, real x0, real gamma) 
	{ 
		return  (1 / (PI * gamma)) * (pow(gamma, 2) / (pow((x - x0), 2) + pow(gamma, 2))); 
	}

	real CauchyCDF(real x, real x0, real gamma) 
	{ 
		return (0.5 + (atan(x - x0) / gamma * PI)); 
	}

	real chiSquareCDF(real k, real x) 
	{ 
		return regularizedGamma(0.5 * k, 0.5 * x); 
	}

	real ChiSquareTest(Matrix A) 
	{ 
		return chiSquareCDF(A.ChiSquareDegreesFreedom(), A.ChiSquareTestStatistic()); 
	}

	real lognormalDistributionPDF(real x, real mean, real variance) 
	{//for non-lognormally distributed data
		//generate log-normal mean and variance
		real mu = log(mean / sqrt(1 + (variance / (mean * mean))));
		real vari = log(1 + (variance / (mean * mean)));

		return (1.0 / (x * sqrt(vari * 2 * PI))) * exp((-0.5 * pow(log(x - mu), 2)) / vari);
	}

	real lognormalDistributionCDF(real x, real mean, real variance) 
	{//for non-lognormally distributed data
		//generate log-normal mean and variance
		real mu = log(mean / sqrt(1 + (variance / (mean * mean))));
		real vari = log(1 + (variance / (mean * mean)));

		return 0.5 * (1 + std::erf(log(x - mu) / sqrt(2 * vari)));
	}

	real FDistributionApproximation(real x, real v1, real v2) 
	{//x is test statistic, v1/v2 are degrees freedom
		real part1 = pow((((((2 * v2) + ((v1 * x) / 3) + v1 - 2) * x) / ((2 * v2) + 
			((4 * v1 * x) / 3)))), 0.333333333333333333333333333333333333333);
		real part2 = sqrt(2 / (9 * v1));
		real part3 = (1 - (2 / (9 * v1)));
		return normalCDF((part1 - part3) / part2);
	}

	real correlationCoefficient(Vector a, Vector b) {
		if (a.size() != b.size()) { return 0; }
		real n = a.size();
		real mu1 = a.mean();
		real mu2 = b.mean();
		real s1 = a.sampleStandardDeviation();
		real s2 = b.sampleStandardDeviation();
		real num = 0;
		real denom1 = 0;
		real denom2 = 0;
		for (int i = 0; i < n; ++i) {
			num += (a[i] - mu1) * (b[i] - mu2);
			denom1 += pow((a[i] - mu1), 2);
			denom2 += pow((b[i] - mu2), 2);
		}
		denom1 = sqrt(denom1);
		denom2 = sqrt(denom2);
		return num / (denom1 * denom2);
	}

	real coefficientOfDetermination(Vector x, Vector y) {
		return pow(correlationCoefficient(x, y), 2);
	}

	Matrix correlationMatrix(Matrix A) {
		//this correlation matrix is the symmetric matrix of values 
		//for Pearson correlations between variable i and variable j
		Matrix M(A.columns, A.columns);
		M.identity();
		for (int i = 1; i < M.columns; ++i) {
			for (int j = 0; j < i; ++j) {
				if (i != j) {
					real temp = correlationCoefficient(A.column(i), A.column(j));
					M(i, j) = temp;
					M(j, i) = temp;
				}
			}
		}
		return M;
	}

	Matrix covarianceMatrix(Matrix M) 
	{//returns a matrix of sample covariances between column vectors
		Matrix A(M.columns, M.columns);
		for (int i = 0; i < M.columns; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				A(i, j) = covariance(M.column(i), M.column(j));
			}
		}
		return A;
	}

	Vector JarqueBeraTest(Matrix M, real criticalValue) 
	{//tests normality
		Vector answer;
		int n = M.rows;
		//real Chi2 = chiSquareCDFInverseApproximation(criticalValue, 2);//Get X^2 value for critical region bound

		for (int i = 0; i < M.columns; ++i) {
			real skewness = 0;
			real kurtosis = 0;
			Vector x = M.column(i);
			real mean = x.mean();
			real s = x.sampleStandardDeviation();

			for (int j = 0; j < n; ++j) {
				skewness += pow(x[j] - mean, 3) / (n * pow(s, 3));
				kurtosis += pow(x[j] - mean, 4) / (n * pow(s, 4));
			}
			kurtosis -= 3;
			real JBtestStat = (n * ((pow(skewness, 2) / 6) + (pow(kurtosis, 2) / 24)));
			answer.push_back(chiSquareCDF(2, JBtestStat));
		}
		return answer;
	}

	Vector JarqueBeraSkewness(Matrix M) 
	{//returns JB skewness vector
		Vector answer;
		int n = M.rows;

		for (int i = 0; i < M.columns; ++i) {
			real skewness = 0;
			Vector x = M.column(i);
			real mean = x.mean();
			real s = x.sampleStandardDeviation();

			for (int j = 0; j < n; ++j) { skewness += pow(x[j] - mean, 3) / (n * pow(s, 3)); }
			answer.push_back(skewness);
		}
		return answer;
	}

	Vector JarqueBeraKurtosis(Matrix M) 
	{//returns JB kurtosis vector
		Vector answer;
		int n = M.rows;

		for (int i = 0; i < M.columns; ++i) {
			real kurtosis = 0;
			Vector x = M.column(i);
			real mean = x.mean();
			real s = x.sampleStandardDeviation();

			for (int j = 0; j < n; ++j) {
				kurtosis += pow(x[j] - mean, 4) / (n * pow(s, 4));
			}
			kurtosis -= 3;
			answer.push_back(kurtosis);
		}
		return answer;
	}

	Vector DAgostinoTest(Matrix M, real criticalValue) 
	{//further refines JarqueBeraTest of normality
		Vector answer;
		real n = M.rows;

		Vector sk = JarqueBeraSkewness(M);
		Vector kt = JarqueBeraKurtosis(M);
		real skewDiv = ((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)));
		real kurtDiv = ((6 * n) / (((n - 2) * (n - 3) * (n + 3) * (n + 5))));
		skewDiv = sqrt(skewDiv);
		kurtDiv = sqrt(kurtDiv) * 2 * (n - 1);

		for (int i = 0; i < sk.size(); ++i) {
			real JBTestStat = pow(sk[i] / skewDiv, 2) + pow(kt[i] / kurtDiv, 2);
			answer.push_back(chiSquareCDF(2, JBTestStat));
		}
		return answer;
	}

	Vector DAgostinoSkewness(Matrix M) 
	{//returns D'Agostino skewness vector
		Vector answer;
		real n = M.rows;

		Vector sk = JarqueBeraSkewness(M);
		real skewDiv = ((6 * n * (n - 1)) / ((n - 2) * (n + 1) * (n + 3)));
		skewDiv = sqrt(skewDiv);

		for (int i = 0; i < sk.size(); ++i) {
			real JBTestStat = sk[i] / skewDiv;
			answer.push_back(JBTestStat);
		}
		return answer;
	}

	Vector DAgostinoKurtosis(Matrix M) 
	{//returns D'Agostino kurtosis vector
		Vector answer;
		real n = M.rows;

		Vector kt = JarqueBeraKurtosis(M);
		real kurtDiv = ((6 * n) / (((n - 2) * (n - 3) * (n + 3) * (n + 5))));
		kurtDiv = sqrt(kurtDiv) * 2 * (n - 1);

		for (int i = 0; i < kt.size(); ++i) {
			real JBTestStat = kt[i] / kurtDiv;
			answer.push_back(JBTestStat);
		}
		return answer;
	}

	Vector AndersonDarlingTest(Matrix M) {
		Vector answer;
		real n = M.rows;

		for (int i = 0; i < M.columns; ++i) {
			Vector x = M.column(i);//data vector
			std::sort(x.begin(), x.end(), [](real a, real b) {return a > b; });
			real mean = x.mean();
			real s = x.sampleStandardDeviation();

			real temp = 0;
			for (int j = 1; j <= M.rows; ++j) {//calculate S
				real z = (x[j - 1] - mean) / s;
				real Yi = normalCDFApproximation(z);
				real part1 = (2 * j) - 1;
				real part2 = (2 * (n - j)) + 1;
				temp += (part1 * log(Yi)) + (part2 * log(1 - Yi));
			}

			//calculate Anderson-Darling test statistic (AD)
			real AD = -n - (temp / n);
			AD *= (1.0 + (0.75 / n) + (2.25 / pow(n, 2)));//adjustment for unknown pop/variance

			//calculate p-value from AD statistic
			real p = 0;
			if (AD >= 0.6) { p = exp(1.2937 - (5.709 * AD) + (0.0186 * pow(AD, 2))); }
			if (0.34 < AD < 0.6) { p = exp(0.9177 - (4.279 * AD) - (1.38 * AD * AD)); }
			if (0.2 < AD <= 0.34) { p = 1 - exp(-8.318 + (42.796 * AD) - (59.938 * AD * AD)); }
			if (AD <= 0.2) { p = 1 - exp(-13.436 + (101.14 * AD) - (223.73 * AD * AD)); }

			//answer.push_back(p);
			answer.push_back(AD);
			x.clear();
		}
		return answer;
	}

	Vector AndersonDarlingTest(Vector X) 
	{
		Vector answer;

		Vector x = X;//data vector
		std::sort(x.begin(), x.end(), [](real a, real b) {return a > b; });
		real mean = x.mean();
		real s = x.sampleStandardDeviation();
		real n = x.size();

		real temp = 0;
		for (int j = 1; j <= x.size(); ++j) {//calculate S
			real z = (x[j - 1] - mean) / s;
			real Yi = normalCDFApproximation(z);
			real part1 = (2 * j) - 1;
			real part2 = (2 * (n - j)) + 1;
			temp += (part1 * log(Yi)) + (part2 * log(1 - Yi));
		}

		// Calculate Anderson-Darling test statistic (AD)
		real AD = -n - (temp / n);
		AD *= (1.0 + (0.75 / n) + (2.25 / pow(n, 2)));//adjustment for unknown pop/variance

		// Calculate p-value from AD statistic
		real p = 0;
		if (AD >= 0.6) { p = exp(1.2937 - (5.709 * AD) + (0.0186 * pow(AD, 2))); }
		if (0.34 < AD < 0.6) { p = exp(0.9177 - (4.279 * AD) - (1.38 * AD * AD)); }
		if (0.2 < AD <= 0.34) { p = 1 - exp(-8.318 + (42.796 * AD) - (59.938 * AD * AD)); }
		if (AD <= 0.2) { p = 1 - exp(-13.436 + (101.14 * AD) - (223.73 * AD * AD)); }

		//answer.push_back(p);
		answer.push_back(AD);
		x.clear();
		return answer;
	}

	real RMS(int T1, int T2, Function f) 
	{
		//format rhs string of f(t) so that it is now abs(f(t)^2)
		std::string funct = "(abs(";
		funct.append(f.function);
		funct.append(")^2)");
		f.function = funct;
		return sqrt(integrateGaussLegendreQuadrature(T1, T2, f) / (T2 - T1));
	}
}
