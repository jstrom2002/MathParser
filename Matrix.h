/**
*   Generic matrix class.
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
		Matrix(std::vector<Vector> elm);
		template <class T>
		Matrix(int row, int col, std::vector<T> element2);	

		bool operator==(Matrix& rhs) { 
			return (element == rhs.element && rows == rhs.rows &&
				columns == rhs.columns);
		}

		//function prototypes
		Matrix add(Matrix A, Matrix B);
		Matrix add(Matrix A, real B);
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
		real maxValue();
		bool canSwapOutZeroDiagonals();
		Matrix characteristicMatrix();
		Polynomial characteristicPolynomial();
		real ChiSquareDegreesFreedom();
		real ChiSquareTestStatistic();
		Matrix Cholesky();
		void clear() { element.clear(); rows = columns = 0; }
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
		Matrix dominantEigenvector();
		real dominantEigenvalue();
		Vector eigenvaluesByGaussianElimination();
		Vector eigenvaluesRealExact();
		Vector eigenvaluesNumerical();
		std::vector<ComplexNumber> eigenvaluesRealAndComplexExact();
		Matrix eigenvalueMatrix();
		Matrix eigenvalueMatrix(Matrix A, int loops);
		Matrix eigenvectors();
		Matrix expandToUpperLeft(int r, int c);
		Matrix expandToLowerRight(int r, int c);
		Matrix extendRows(int r);
		Matrix extendColumns(int c);
		Matrix exponent(Matrix A, int b);
		Vector findRow(Vector x);
		Vector findColumn(Vector x);
		real FTestStatistic();
		real FrobeniusNorm();
		Matrix GaussianElimination();//exact calculation of row-echelon form
		real geometricMean();
		real get(int i, int j);
		int getPivot(int rw);
		int getReversePivot(int rw);
		Matrix GivensRotationMatrix(real a1, real a2, int r1, int r2);
		Matrix GramMatrix();
		Matrix GramMatrix(Matrix M);
		Matrix GramSchmidt();
		Matrix HadamardProduct(Matrix A, Matrix B);
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
		bool isOverDetermined();
		bool isInconsistent(Vector b);
		bool isPositiveDefinite();
		bool isPositiveSemidefinite();
		bool isSingular();
		bool isSymmetric();
		bool isUnderDetermined();
		std::vector<int> JacobiIndexing();
		Matrix JacobiRotationMatrix(int p, int q, real c, real s);
		std::vector<Matrix> JacobiTransformation();
		Matrix L();
		std::vector<Matrix> LDL();
		int largestElement();
		Matrix leastSquares(Matrix X, Matrix Y);
		Function leastSquares();
		Matrix leastSquaresMatrix();
		Matrix lowerTriangularize();
		Matrix lowerHessenbergForm();
		std::vector<Matrix> LU();
		real mean();
		Vector meanVector();
		Matrix multiply(Matrix A, Matrix B);
		Matrix multiply(Matrix A, real B);
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
		Matrix operator*=(Matrix rhs);
		Matrix operator*=(real x);
		Matrix operator*(Matrix rhs);
		Matrix operator*(real x);
		Matrix operator+=(Matrix rhs);
		Matrix operator+=(real x);
		Matrix operator+(Matrix rhs);
		Matrix operator+(real x);
		Matrix operator-=(Matrix rhs);
		Matrix operator-=(real x);
		Matrix operator-(Matrix rhs);
		Matrix operator-(real x);
		bool operator==(Matrix B);
		Matrix outerProduct(Matrix A, Matrix B);
		real pNorm(real p);
		std::vector<int> pivotColumns();
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
		void set(int i, int j, real x);
		void setColumnNonZeroValues(int col, real x);
		void setRowNonZeroValues(int rw, real x);
		void setRow(int rw, Vector n);
		void setColumn(int col, Vector n);
		int size();
		bool solutionCheck(Vector x, Vector b);
		Vector solve(Vector b);
		Matrix submatrix(int i, int j, int i2, int j2);
		Matrix subtract(Matrix A, Matrix B);
		Matrix subtract(Matrix A, real B);
		real sum();
		real sumAll();
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
		std::vector<Vector> toVectorArray();
		Matrix transpose();
		real trace();
		void trim();
		Matrix U();
		Matrix upperTriangularize();
		Matrix upperHessenbergForm();
		Matrix Vandermonde();
		Matrix varianceMatrix();
		bool zeroInDiag();
	};

	Matrix directionMatrix(Vector v);
	Matrix identityMatrix(int n);
	Matrix identityMatrix(int n, int m);
	Matrix positionMatrix(Vector v);
	Matrix rotationMatrix(Vector v);
	Matrix scalingMatrix(Vector v);
	Matrix translationMatrix(Vector v);
	Matrix Vandermonde(Polynomial p);
}

namespace std
{
	MathParser::Matrix pow(MathParser::Matrix m, MathParser::real x);
}