/**
*   Generic matrix class. All matrices contain real valued members. Entries are stored in
*	an array indexed by row-major indexing. Accessing entries is performed by use of the
*	'()' operator, ie 'A(0,0) = 1.'
*/

#pragma once
#include "Types.h"
#include "BaseObject.h"
#include "Vector.h"
#include <vector>

namespace MathParser
{
	class ComplexNumber;
	class Function;
	class Polynomial;
	class Vector;

	class Matrix : public BaseObject
	{
	public:
		Vector element;
		int rows;
		int columns;

		Matrix();
		Matrix(int n);
		Matrix(int row, int col);
		Matrix(int row, int col, Vector elements);
		Matrix(Vector elements);
		Matrix(std::vector<Vector> elements);
		Matrix(int row, int col, std::vector<real> elements);	

		// Overloaded operators:
		bool operator==(Matrix& rhs);
		bool operator!=(Matrix& rhs);
		real& operator()(int r, int c);
		real operator()(int r, int c) const;
		Matrix operator*=(Matrix& rhs);
		Matrix operator*=(real x);
		Matrix operator+=(Matrix& rhs);
		Matrix operator+=(real x);
		Matrix operator+(Matrix& rhs);
		Matrix operator+(real x);
		Matrix operator-=(Matrix& rhs);
		Matrix operator-=(real x);
		friend Matrix operator+(Matrix& lhs, Matrix& rhs);
		friend Matrix operator+(Matrix& lhs, real x);
		friend Matrix operator+(real x, Matrix& rhs);
		friend Matrix operator-(Matrix& lhs, Matrix& rhs);
		friend Matrix operator-(Matrix& lhs, real x);
		friend Matrix operator-(real x, Matrix& rhs);
		friend Matrix operator*(Matrix& lhs, Matrix& rhs);
		friend Matrix operator*(Matrix& lhs, real x);
		friend Matrix operator*(real x, Matrix& rhs);
		friend Matrix operator/(Matrix& lhs, Matrix& rhs);
		friend Matrix operator/(Matrix& lhs, real x);
		friend Matrix operator/(real x, Matrix& rhs);

		// Class method prototypes:
		Matrix addColumn(Vector vec);
		Matrix addColumn(int j, Vector vec);
		Matrix addRow(Vector vec);
		Matrix addRow(int i, Vector vec);
		void addToRow(int rw, real n);
		void addToRow(int rw, Vector n);
		void addToColumn(int cl, real n);
		void addToColumn(int cl, Vector n);
		real adjustedRsquared();
		Matrix characteristicMatrix();
		Polynomial characteristicPolynomial();
		real ChiSquareDegreesFreedom();
		real ChiSquareTestStatistic();
		Matrix Cholesky();
		void clear();
		Vector column(int c);
		real columnAbsMax(int r);
		real columnMax(int r);
		real columnCovariance(int c1, int c2);
		real columnDeviation(int c);
		real columnMean(int c);
		real columnAbsMin(int r);
		real columnMin(int r);
		int columnNonzeroValues(int rw);
		real columnNorm(int c);
		real columnSumSquares(int c);
		real columnVariance(int c);
		real det();
		Vector diagonal();
		Matrix directSum(Matrix A, Matrix B);
		Matrix dominantEigenvector(int iterations = 500);
		real dominantEigenvalue(int iterations = 500);
		Vector eigenvaluesByGaussianElimination();
		Vector eigenvaluesRealExact();
		Vector eigenvaluesNumerical();
		std::vector<ComplexNumber> eigenvaluesRealAndComplexExact(int iterations=10000);
		Matrix eigenvalueMatrix();
		Matrix eigenvalueMatrix(Matrix A, int loops);
		Matrix eigenvectors(int iterations = 500);
		Matrix expandToUpperLeft(int r, int c);
		Matrix expandToLowerRight(int r, int c);
		Matrix extendRows(int r);
		Matrix extendColumns(int c);
		Vector findRow(Vector x);
		Vector findColumn(Vector x);
		real FTestStatistic();
		real FrobeniusNorm();
		Matrix GaussianElimination();//exact calculation of row-echelon form
		real geometricMean();
		int getPivot(int rw);
		int getReversePivot(int rw);
		Matrix GivensRotationMatrix(real a1, real a2, int r1, int r2);
		Matrix GramMatrix();
		Matrix GramMatrix(Matrix M);
		Matrix GramSchmidt();
		std::vector<Matrix> HouseholderQR();
		Matrix inverse();
		Matrix inverseExact();
		Matrix inverseByQR();
		bool isBoolean();
		bool isDiagonalized();
		bool isIntegerMatrix();
		bool isLinearlyIndependent();
		bool isOrthonormal();
		bool isInconsistent(Vector b);
		bool isPositiveDefinite();
		bool isPositiveSemidefinite();
		bool isSymmetric();
		std::vector<Matrix> LDL();
		Matrix leastSquares(Matrix X, Matrix Y);
		Function leastSquares();
		Matrix leastSquaresMatrix();
		Matrix lowerTriangularize();
		Matrix lowerHessenbergForm();
		real max();
		real mean();
		Vector meanVector();
		real min();
		void multiplyRow(int rw, real n);
		void multiplyRow(int rw, Vector n);
		real norm();
		void normalize();
		void normalizeColumn(int c);
		void normalizeRow(int r);
		Vector normalizedColumn(int c);
		Vector normalizedRow(int r);
		int nullity();
		Matrix nullSpace();
		real pNorm(real p);
		Matrix populationCovarianceMatrix();
		real populationStandardDeviation();
		Matrix pseudoinverse();
		std::vector<Matrix> QR();
		void randomize();
		void randomizeInteger();
		void randomizeBoolean();
		void randomizeSymmetric();
		Matrix reducedRowEchelonForm();//numerical method of calculating RRE form
		void removeColumn(int n);
		void removeRow(int n);
		void removeZeroRows();
		void removeZeroColumns();
		void reshape(int r, int c);
		Vector residuals();
		void reverseColumn(int n);
		void reverseRow(int n);
		Vector row(int r);
		real rowCovariance(int r1, int r2);
		real rowDeviation(int r);
		real rowNorm(int r);
		real rowAbsMax(int r);
		real rowMax(int r);
		real rowMean(int r);
		real rowAbsMin(int r);
		real rowMin(int r);
		int rowNonzeroValues(int rw);
		int rank();
		real rowSumSquares(int r);
		real rowVariance(int r);
		real Rsquared();
		Matrix sampleCovarianceMatrix();
		real sampleStandardDeviation();
		void setRow(int rw, Vector n);
		void setColumn(int col, Vector n);
		int size();
		Vector solve(Vector b);
		Matrix submatrix(int i, int j, int i2, int j2);
		real sum();
		real sumColumn(int c);
		real sumRow(int r);
		real sumSquares();
		void swapRows(int r1, int r2);
		void swapColumns(int c1, int c2);
		std::vector<Matrix> SVD();
		Matrix tensorProduct(Matrix A, Matrix B);
		std::vector<Matrix> thinQR();
		std::vector<Matrix> thinHouseholderQR();
		virtual std::string to_string(int precision = 4);
		Matrix transpose();
		real trace();
		void trim();
		Matrix U();
		Matrix upperTriangularize();
		Matrix upperHessenbergForm();
		Matrix Vandermonde();
		Matrix varianceMatrix();

		private:
			// Helper functions to reduce reusing code for operator overloads.
			Matrix add(Matrix& a, Matrix& b);
			Matrix add(Matrix& a, real  b);
			Matrix add(real b, Matrix& a);
			Matrix subtract(Matrix& a, Matrix& b);
			Matrix subtract(Matrix& a, real b);
			Matrix subtract(real b, Matrix& a);
			Matrix multiply(Matrix& a, Matrix& b);
			Matrix multiply(Matrix& a, real b);
			Matrix multiply(real b, Matrix& a);
			Matrix divide(Matrix& p, Matrix& q);
			Matrix divide(Matrix& p2, real q);
			Matrix divide(real q, Matrix& p2);		
	};

	// Non-class member functions for matrices:
	Matrix hconcat(Matrix A, Matrix B);
	Matrix vconcat(Matrix A, Matrix B);
	Matrix convolve(Matrix& kernel, Matrix& B);
	Matrix crossCorrelation(Matrix A, Matrix B);
	Matrix csvread(std::string filename);
	void csvwrite(std::string filename, Matrix& M);
	Matrix directionMatrix(Vector v);
	Matrix HadamardProduct(Matrix A, Matrix B);
	Matrix eye(int n);
	Matrix eye(int n, int m);
	Matrix imread(std::string filename);
	void imwrite(std::string filename, Matrix& M, int colorDepth=8);
	std::vector<Matrix> JacobiTransformation(Matrix& M, highpUint limit = 10000);
	std::vector<Matrix> LU(Matrix& M);
	Matrix outerProduct(Matrix& A, Matrix& B);
	Matrix positionMatrix(Vector v);
	Matrix rotationMatrix(Vector v);
	Matrix scalingMatrix(Vector v);
	Matrix translationMatrix(Vector v);
	Matrix Vandermonde(Polynomial p);
	Matrix wedgeProduct(Matrix& A, Matrix& B);
}

namespace std
{
	MathParser::Matrix pow(MathParser::Matrix& m, int x);
}