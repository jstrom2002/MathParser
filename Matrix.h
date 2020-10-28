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
		Matrix(int row, int col, Vector element2);
		Matrix(Vector element2);
		Matrix(std::vector<Vector> elm);
		Matrix(int row, int col, std::vector<real> element2);	

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
		Matrix autoCorrelation(Matrix A);
		bool canSwapOutZeroDiagonals();
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
		Matrix concatenate(Matrix A, Matrix B);
		Matrix convolve(Matrix A, Matrix B);
		Matrix crossCorrelation(Matrix A, Matrix B);
		real det();
		real detExact();
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
		void identity();
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
		std::vector<int> JacobiIndexing();
		Matrix JacobiRotationMatrix(int p, int q, real c, real s);
		std::vector<Matrix> JacobiTransformation(highpUint limit = 10000);
		Matrix L();
		std::vector<Matrix> LDL();
		int largestElement();
		Matrix leastSquares(Matrix X, Matrix Y);
		Function leastSquares();
		Matrix leastSquaresMatrix();
		void loadBMP(std::string filename);
		void loadCSV(std::string filename);
		Matrix lowerTriangularize();
		Matrix lowerHessenbergForm();
		std::vector<Matrix> LU();
		real maxValue();
		real mean();
		Vector meanVector();
		real minValue();
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
		void power(real n);
		Matrix pseudoinverse();
		Matrix Q();
		std::vector<Matrix> QR();
		Matrix R();
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
		void saveBMP(std::string filename, int colorDepth=8);
		void saveCSV(std::string filename);
		void setColumnNonZeroValues(int col, real x);
		void setRowNonZeroValues(int rw, real x);
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
		Matrix swapOutZeroDiagonals();
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
		bool zeroInDiag();

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
			
			// Helper functions for class methods:
			std::vector<int> pivotColumns();
	};

	Matrix convolve(Matrix kernel, Matrix img);
	Matrix directionMatrix(Vector v);
	Matrix identityMatrix(int n);
	Matrix identityMatrix(int n, int m);
	Matrix HadamardProduct(Matrix A, Matrix B);
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
	MathParser::Matrix pow(MathParser::Matrix m, int x);
}