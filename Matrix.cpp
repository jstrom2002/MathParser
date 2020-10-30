#include "Matrix.h"
#include "Polynomial.h"
#include "Function.h"
#include "Vector.h"
#include "ComplexNumber.h"
#include "MathLib.h"
#include "StringUtils.h"
#include "Parsing.h"
#include "FileIO.h"
#include <ctime>
#include <iomanip>
#include <fstream>

namespace std
{
	MathParser::Matrix pow(MathParser::Matrix& A, int b)
	{
		if (b == 0)
			return MathParser::eye(A.rows, A.columns);

		MathParser::Matrix A2 = A;
		if (b < 0)
		{
			A2 = A.inverse();
			b = abs(b);
		}

		for (int i = 1; i < b; ++i)		
			A2 *= A;		

		return A2;
	}
}

namespace MathParser
{
	Matrix::Matrix() : rows(0), columns(0)
	{ 
	}

	Matrix::Matrix(int n) 
	{
		rows = n;
		columns = n;
		element.resize(n * n, 0);
	}

	Matrix::Matrix(int row, int col) 
	{
		element.resize(row * col, 0);
		rows = row; 
		columns = col;
	}

	Matrix::Matrix(int row, int col, Vector elements) 
	{
		int n = row;
		int m = col;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row * col, 0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				this->element[i * col + j] = elements[i * col + j];
			}
		}
		rows = row; 
		columns = col;
	}

	Matrix::Matrix(Vector v)
	{
		rows = columns = v.size();
		element.resize(v.size() * v.size(), 0);
		for (int i = 0; i < v.size(); ++i)
			element[i * v.size() + i] = v[i];
	}

	Matrix::Matrix(int row, int col, std::vector<real> elements) 
	{
		int n = row;
		int m = col;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row * col, 0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				element[i * col + j] = elements[i * col + j];				
			}
		}
		rows = row; 
		columns = col;
	}

	Matrix::Matrix(std::vector<Vector> elements) 
	{
		rows = elements.size();
		columns = elements[0].size();
		element.clear();
		element.resize(rows*columns, 0);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				element[i * columns + j] = elements[i][j];
			}
		}
	}

	bool zeroInDiag(Matrix& M)
	{	//Test if there are any 0's in the diagonal of a matrix.
		for (int i = 0; i < M.rows; ++i)
			if (!M(i, i))
				return true;
		return false;
	}

	std::vector<int> pivotColumns(Matrix& M) 
	{	// Helper function.
		std::vector<int> piv;
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				real val = M(i, j);
				if (val != 0) {
					piv.push_back(j);
					j = M.columns + 1;
				}
			}
		}
		return piv;
	}

	void setRowNonZeroValues(Matrix& M, int rw, real x)
	{	// Helper function.
		for (int k = 0; k < M.columns; ++k) {
			if (std::abs(M(rw, k)) > 0) {
				M(rw, k) = x;
			}
		}
	}

	void setColumnNonZeroValues(Matrix& M, int col, real x)
	{	// Helper function.
		for (int k = 0; k < M.rows; ++k)
		{
			if (std::abs(M(k, col)) > 0)
				M(k, col) = x;
		}
	}

	real& Matrix::operator()(int r, int c) 
	{ 
		return element[r * columns + c]; 
	}

	real Matrix::operator()(int r, int c) const 
	{ 
		return element[r * columns + c]; 
	};

	void Matrix::clear() 
	{ 
		element.clear(); 
		rows = columns = 0; 
	}

	int Matrix::size() 
	{ 
		return rows * columns; 
	}

	real Matrix::max() 
	{
		real x = -std::numeric_limits<real>::max();
		for (int i = 0; i < rows; ++i) 
		{
			for (int j = 0; j < columns; ++j) 
			{
				if (element[i * columns + j] > x) 
					x = element[i * columns + j];
			}
		}
		return x;
	}

	real Matrix::min()
	{
		real x = std::numeric_limits<real>::max();
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < columns; ++j)
			{
				if (element[i * columns + j] < x)
					x = element[i * columns + j];
			}
		}
		return x;
	}

	Matrix Matrix::submatrix(int i, int j, int i2, int j2) 
	{
		int size1 = (std::abs(i2-i)+1) * (std::abs(j2-j)+1);
		Vector vals;
		vals.resize(size1,0);
		int counter = 0;
		for (int a = i; a <= i2; ++a) 
			for (int b = j; b <= j2; ++b) {
				vals[counter] = element[a * columns + b];
				++counter;
			}		
		return Matrix(std::abs(i2 - i)+1, std::abs(j2 - j)+1, vals);
	}

	Vector Matrix::diagonal() 
	{
		Vector vals;
		for (int i = 0; i < rows; ++i) {
			vals.push_back(element[i * columns + i]);
		}
		return vals;
	}

	Matrix Matrix::inverseExact() 
	{
		if (rows != columns)
			return Matrix();
		Matrix m(rows, columns, element);
		m = hconcat(m, eye(rows, columns));
		m = m.GaussianElimination();
		while (m.columns > columns)
			m.removeColumn(0);
		return m;
	}

	Matrix Matrix::inverse() {
		/*INPUT: n×m matrix A.
		OUTPUT: n×m matrix in reduced row echelon form.
		1. Set j <-- 1
		2. For each row i from 1 to n do
		a. While column j h as all zero elements, set j <- j+1. If j>m return A.
		b. If element a_ij is zero, then interchange row i with a row x > i that has a_xj!=0.
		c. Divide each element of row i by a_ij, thus making the pivot aij equal to one.
		d. For each row k from 1 to n, with k != i, subtract row i multiplied by a_kj from row k.
		3. Return transformed matrix A.		*/

		if (det() == 0) 
			return Matrix();
		Matrix A(rows, columns, element);
		Matrix Aug = eye(rows, columns);

		while (zeroInDiag(A)) {
			for (int i = 0; i < rows; ++i) {
				if (A(i, i) == 0) {
					for (int j = 0; j < columns; ++j) {
						if (A(i, j) != 0 && A(j, i) != 0 && i != j) 
						{ //choose rows to swap to create a non-zero diagonal
							A.swapRows(i, j);
							Aug.swapRows(i, j);
							j = columns;
						}
					}
				}
			}
		}

		//go down and eliminate the lower triangle of values
		for (int j = 0; j < columns; ++j) 
		{
			real pivot = A(j, j);
			for (int i = 0; i < rows; ++i) {
				if (i != j && A(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A(i, j)) == abs(pivot)) {
						if (A(i, j) == pivot) 
						{
							temp = temp * -1.0;
							temp2 = temp2 * -1.0;
							goto skip;
						}
						if ((A(i, j) == (-1 * pivot)) || ((-1 * A(i, j)) == pivot)) 
							goto skip;
					}
					if (abs(A(i, j)) != abs(pivot)) {
						real a = -A(i, j) / pivot;
						temp = temp * a;
						temp2 = temp2 * a;
					}

				skip:
					A.addToRow(i, temp);
					Aug.addToRow(i, temp2);

				}
			}
		}

		//go the other way to remove the upper triangle of non-diagonal values
		for (int j = columns - 1; j >= 0; --j) 
		{
			real pivot = A(j, j);
			for (int i = rows - 1; i >= 0; --i) {
				if (i != j && A(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A(i, j)) == abs(pivot)) {
						if (A(i, j) == pivot) {
							temp = temp * -1.0;
							temp2 = temp2 * -1.0;
							goto skip2;
						}
						if ((A(i, j) == (-1 * pivot)) || ((-1 * A(i, j)) == pivot)) { goto skip2; }
					}
					if (abs(A(i, j)) != abs(pivot)) {
						real a = -A(i, j) / pivot;
						temp = temp * a;
						temp2 = temp2 * a;
					}

				skip2:
					A.addToRow(i, temp);
					Aug.addToRow(i, temp2);
				}
			}
		}

		for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
			if (A(i, i) != 0 && A(i, i) != 1) {
				real a = ((real)1 / (A(i, i)));
				A.multiplyRow(i, a);
				Aug.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, Aug.element);
	}

	Matrix HadamardProduct(Matrix A, Matrix B) 
	{// Multiplies elements of two matrices by index so that C(i,j) = A(i,j)*B(i,j)
		if (A.rows != B.rows || A.columns != B.columns) 
			return Matrix();
		Vector elm(A.size());
		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < A.columns; ++j) 
			{
				elm[i * A.columns + j] = A(i, j) * B(i, j);
			}
		}
		return Matrix(A.rows, A.columns, elm);
	}

	void Matrix::reshape(int r, int c)
	{
		// Catch case where new dimensions are invalid.
		if (r * c != rows * columns)
			return;

		else
		{
			rows = r;
			columns = c;
		}
	}

	void imwrite(std::string filename, Matrix& M, int colorDepth)
	{//Ref: https://en.wikipedia.org/wiki/BMP_file_format#Example_1
		std::string ext = getExtension(filename);
		if (ext == ".bmp")
			saveBMP(filename, M.rows, M.columns, 24, M.element.get());
		else if (ext == ".ppm")
			savePPM(filename, M.rows, M.columns, 24, M.element.get());
		else if (ext == ".tga")
			saveTGA(filename, M.rows, M.columns, 24, M.element.get());
	}

	Matrix imread(std::string filename)
	{
		std::vector<real> output;
		int rows = 0;
		int columns = 0;
		std::string ext = getExtension(filename);		

		if (ext == ".bmp")
			loadBMP(filename, &rows, &columns, &output);
		else if (ext == ".ppm")
			loadPPM(filename, &rows, &columns, &output);
		else if (ext == ".tga")
			loadTGA(filename, &rows, &columns, &output);

		// Add empty pixels if there was some error during image loading.
		while (output.size() < rows * columns)
			output.push_back(0);

		return Matrix(rows, columns, output);
	}

	Matrix convolve(Matrix kernel, Matrix img)
	{
		// Handle case where kernel has an even dimension.
		if (!(kernel.rows % 2) || !(kernel.columns % 2))
			return Matrix();

		int dx = std::floor(kernel.columns/2);
		int dy = std::floor(kernel.rows/2);
		std::vector<real> vals(img.size(), 0);
		for(int i = dy; i < img.rows - dy; ++i)
			for (int j = dx; j < img.columns - dx; ++j)
			{// Each region is convolved by taking the Hadamard product of it with the kernel
			 // and summing the results.
				Matrix submat = img.submatrix(i - dy, j - dx, i + dy, j + dx);
				Matrix had = HadamardProduct(submat, kernel);
				vals[i * img.columns + j] = had.sum();
			}
		return Matrix(img.rows, img.columns, vals);
	}

	Vector Matrix::row(int rw) 
	{
		if (rw > rows || rw < 0) 		
			return Vector(); 
		
		Vector x;
		x.resize(columns,0);
		for (int k = 0; k < columns; ++k)
		{
			x[k] = element[rw * columns + k];
		}
		return x;
	}

	Vector Matrix::column(int cl) 
	{
		if (cl > columns || cl < 0)
			return Vector();

		Vector x;
		x.resize(rows,0);
		for (int k = 0; k < rows; ++k) {
			x[k] = element[k * columns + cl];
		}
		return x;
	}

	void Matrix::setRow(int rw, Vector n) 
	{
		for (int k = 0; k < columns; ++k) 
			element[rw * columns + k] = n[k];
	}

	void Matrix::setColumn(int col, Vector n) {
		for (int k = 0; k < rows; ++k) 
		{ 
			element[k * columns + col] = n[k];
		}
	}

	void Matrix::multiplyRow(int rw, real n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(rw * columns) + k] *= n;
		}
	}

	void Matrix::multiplyRow(int rw, Vector n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(rw * columns) + k] *= n[k];
		}
	}

	void Matrix::addToRow(int rw, real n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(rw * columns) + k] += n;
		}
	}

	void Matrix::addToRow(int rw, Vector n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(rw * columns) + k] += n[k];
		}
	}

	void Matrix::addToColumn(int cl, real n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(k * columns) + cl] += n;
		}
	}

	void Matrix::addToColumn(int cl, Vector n) {
		for (int k = 0; k < columns; ++k)
		{
			element[(k * columns) + cl] += n[k];
		}
	}

	Vector Matrix::findRow(Vector x) 
	{
		int index = 0;
		bool isRow = false;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (x[j] != element[i * columns + j]) 
				{ 
					isRow = false; 
					j = columns; 
				}
				if (j = columns - 1) 
				{ 
					isRow = true; 
					index = i; 
				}
			}
		}
		return row(index);
	}

	Vector Matrix::findColumn(Vector x) 
	{
		int index = -1;
		bool isColumn = false;
		for (int j = 0; j < columns; ++j) {
			for (int i = 0; i < rows; ++i) {
				if (x[j] != element[i * columns + j]) 
				{ 
					isColumn = false; 
					i = rows; 
				}
				if (i = rows - 1) 
				{ 
					isColumn = true; 
					index = j; 
				}
			}
		}
		if (index > -1) 
			return column(index);
		else 
		{ 
			Vector n; 
			return n; 
		}
	}

	void Matrix::reverseRow(int n) 
	{
		Vector vec = row(n);
		std::reverse(vec.begin(), vec.end());
		for (int j = 0; j < columns; ++j) {
			element[n * columns + j] = vec[j];
		}
	}

	void Matrix::reverseColumn(int n) 
	{
		Vector vec = column(n);
		std::reverse(vec.begin(), vec.end());
		for (int j = 0; j < rows; ++j) {
			element[j * columns + n] = vec[j];
		}
	}

	void Matrix::swapRows(int rw1, int rw2) {
		Vector r1 = row(rw1);
		Vector r2 = row(rw2);
		for (int i = 0; i < rows; ++i) 
		{
			for (int j = 0; j < columns; ++j) 
			{
				if (i == rw1) 
					element[i * columns + j] = r2[j];
				if (i == rw2) 
					element[i * columns + j] = r1[j];
				if (i > rw1 && i > rw2) 
					return;
			}
		}
	}

	void Matrix::swapColumns(int col1, int col2) 
	{
		Vector c1 = column(col1);
		Vector c2 = column(col2);
		for (int j = 0; j < columns; ++j) {
			for (int i = 0; i < rows; ++i) {
				if (j == col1) 
					element[i * columns + j] = c2[i];
				if (j == col2) 
					element[i * columns + j] = c1[i];
				if (j > col1 && j > col2) 
					return;
			}
		}
	}

	void Matrix::randomize() 
	{
		srand(clock());
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real a = ((real)rand() / (pow(10, (rand() % 5) + 3))) * 
					pow(-1, rand() % 2);
				while (a == 0) { ((real)rand() / (pow(10, (rand() % 5) + 3)))* pow(-1, rand() % 2); }
				element[i * columns + j] = a;
			}
		}
	}

	void Matrix::randomizeInteger() 
	{
		srand(clock());
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real a = (int)((real)rand() / (pow(10, (rand() % 4) + 2))) * 
					pow(-1, rand() % 2);
				element[i * columns + j] = a;
			}
		}
	}

	void Matrix::randomizeBoolean() 
	{
		srand(clock());
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real a = 1 * (rand() % 2);
				element[i * columns + j] = a;
			}
		}
	}

	void Matrix::randomizeSymmetric() 
	{
		srand(clock());
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j <= i; ++j) {
				element[i * columns + j] = ((real)rand() / (30000)) * 
					pow(-1, rand() % 2);
				element[j * columns + i] = element[i * columns + j];
			}
		}
	}

	Matrix Matrix::expandToUpperLeft(int r, int c) 
	{//takes original matrix and creates a new matrix with the original 
	 //in the lower right corner.
		Matrix newEl(r, c);
		if (r > rows && c > columns) {
			newEl = eye(r,c);
			for (int i = r - rows; i < r; ++i) {
				for (int j = c - columns; j < c; ++j) {
					int idx_i = (i - (r - rows));
					int idx_j = (j - (c - columns));
					newEl(i, j) = element[idx_i * columns + idx_j];
				}
			}
		}
		return newEl;
	}

	Matrix Matrix::expandToLowerRight(int r, int c) 
	{//takes original matrix and creates a new matrix with the original 
	 //in the upper left corner.
		Matrix newEl(r, c);
		if (r > rows && c > columns) {
			newEl = eye(r, c);
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < columns; ++j) {
					newEl(i, j) = element[i * columns + j];
				}
			}
		}
		return newEl;
	}

	Matrix Matrix::extendRows(int r) {
		Matrix newEl(rows + r, columns);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				newEl(i, j) = element[i * columns + j];
			}
		}
		return Matrix(rows + r, columns, newEl.element);
	}

	Matrix Matrix::extendColumns(int c) {
		Matrix newEl(rows, columns + c);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				newEl(i, j) = element[i * columns + j];
			}
		}
		return Matrix(rows, columns + c, newEl.element);
	}

	Matrix hconcat(Matrix A, Matrix B)
	{
		int r = A.rows;
		if (B.rows > A.rows) 
			r = B.rows;
		int c = A.columns + B.columns;
		Matrix C(r, c);
		for (int i = 0; i < r; ++i) 
		{
			for (int j = 0; j < c; ++j) 
			{
				if (j < A.columns && i < A.rows) 
					C(i, j) = A(i, j);
				if (j >= A.columns && i < B.rows) 
					C(i, j) = B(i, j - A.columns);
			}
		}
		return C;
	}

	Matrix vconcat(Matrix A, Matrix B)
	{
		int r = A.rows + B.rows;
		int c = A.columns;
		if (B.columns > A.columns)
			c = B.columns;
		Matrix C(r, c);
		for (int i = 0; i < r; ++i)
		{
			for (int j = 0; j < c; ++j)
			{
				if (i < A.rows && j < A.columns)
					C(i, j) = A(i, j);
				if (i >= A.rows && j < B.columns)
					C(i, j) = B(i - A.rows, j);
			}
		}
		return C;
	}

	Matrix crossCorrelation(Matrix A, Matrix B) 
	{
		int kCenterX = A.columns / 2;
		int kCenterY = A.rows / 2;

		int mm, nn, ii, jj;

		Vector vals(B.element.size());

		for (int i = kCenterY; i < B.rows-kCenterY; ++i) {
			for (int j = kCenterX; j < B.columns-kCenterX; ++j) {
				for (int m = 0; m < A.rows; ++m) 
				{
					for (int n = 0; n < A.columns; ++n) 
					{
						// Index of input signal, used for checking boundaries.
						ii = i + (kCenterY - m);
						jj = j + (kCenterX - n);

						// Ignore input samples which are out of bounds.
						if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) 
							vals[i * B.columns + j] += B(ii, jj) * A(m, n);
					}
				}
			}
		}
		return Matrix(B.rows, B.columns, vals);
	}


	Matrix convolve(Matrix& kernel, Matrix& B) 
	{	// Covolve Matrix kernel w/ Matrix B
		int kCenterX = kernel.columns / 2;
		int kCenterY = kernel.rows / 2;

		int mm, nn, ii, jj;

		Vector vals(B.element.size());

		for (int i = kCenterY; i < B.rows - kCenterY; ++i) 
		{
			for (int j = kCenterX; j < B.columns - kCenterX; ++j) 
			{
				for (int m = 0; m < kernel.rows; ++m) 
				{
					mm = kernel.rows - 1 - m;// row index of flipped kernel
					for (int n = 0; n < kernel.columns; ++n) 
					{
						nn = kernel.columns - 1 - n;// column index of flipped kernel

						// index of input signal, used for checking boundary
						ii = i + (kCenterY - mm);
						jj = j + (kCenterX - nn);

						// ignore input samples which are out of bounds.
						if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) 
							vals[i * B.columns + j] += B(ii, jj) * kernel(mm, nn);
					}
				}
			}
		}
		return Matrix(B.rows, B.columns, vals);
	}

	bool Matrix::isIntegerMatrix() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) 
			{
				real val = element[i * columns + j];
				if (std::floor(std::abs(val)) !=  std::abs(val))
					return false;
			}
		}
		return true;
	}

	bool Matrix::isSymmetric() {
		Matrix M = transpose();
		if (*this == M) 
			return true;
		return false;
	}

	bool Matrix::isBoolean() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real val = element[i * columns + j];
				if (val && val != 1) 
					return false;
			}
		}
		return true;
	}

	bool Matrix::isDiagonalized() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i != j) {
					if (std::abs(element[i * columns + j]) > 0)
						return false;
				}
			}
		}
		return true;
	}

	bool Matrix::isOrthonormal() 
	{
		for (int i = 0; i < columns; ++i) 
		{
			if (columnNorm(i) != 1)
				return false;
		}
		if (!det()) 
			return false;		
		return true;
	}

	bool Matrix::isLinearlyIndependent() {
		if (rows < columns) 
			return true;
		if (!det()) 
			return false;
		if (!GramMatrix().det()) 
			return false;
		return true;
	}

	Vector Matrix::meanVector() 
	{
		Vector vec;
		for (int i = 0; i < columns; ++i) {
			vec.push_back(columnMean(i));
		}
		return vec;
	}

	Matrix Matrix::varianceMatrix() 
	{
		Matrix M(columns, columns);
		for (int i = 0; i < columns; ++i) 
		{
			for (int j = 0; j < rows; ++j) 
			{
				element[i*rows + j] = element[j * columns + i];
			}
		}
		return Matrix();
	}

	Matrix Matrix::populationCovarianceMatrix() 
	{
		Matrix M(columns, columns);
		for (int i = 0; i < columns; ++i) 
		{
			for (int j = 0; j < columns; ++j) 
			{
				Vector temp_i = column(i);
				Vector temp_j = column(j);
				Vector X1 =  temp_i - columnMean(i);
				Vector X2 = temp_j - columnMean(j);
				M(i, j) = dot(X1, X2) / (rows);
			}
		}
		return M;
	}

	Matrix Matrix::sampleCovarianceMatrix() 
	{
		Matrix M(columns, columns);
		for (int i = 0; i < columns; ++i) {
			for (int j = 0; j < columns; ++j) {
				Vector temp_i = column(i);
				Vector temp_j = column(j);
				Vector X1 = temp_i - columnMean(i);
				Vector X2 = temp_j - columnMean(j);
				M(i, j) = dot(X1, X2) / ((rows)-1);
			}
		}
		return M;
	}

	Matrix Matrix::transpose() 
	{
		int n = rows;
		int m = columns;
		Vector p(n * m);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				p[i * rows + j] = element[j * columns + i];
			}
		}
		return Matrix(m, n, p);
	}

	Matrix Matrix::GramMatrix() 
	{//turn row vector or nxn matrix into its Gram matrix
		Matrix T = transpose();
		if (columns == 1) 
			return *this * T;
		else 
			return T * *this;
	}

	Matrix Matrix::GramMatrix(Matrix M) 
	{//turn 2 vectors into a Gram matrix
		Matrix Trans = M.transpose();
		if (M.columns == 1) 
			return M * Trans;
		else 
			return Trans * M;
	}

	Matrix Matrix::Vandermonde() 
	{//create Vandermonde matrix from vector
		int p = columns;
		if (columns == 1) 
			p = rows;
		Vector temp;
		temp.resize(p * p,0);
		for (int i = 0; i < p; ++i) {
			for (int j = 0; j < p; ++j) {
				temp[(i * columns) + j] = pow(element[i], j);
			}
		}
		return Matrix(p, p, temp);
	}

	real Matrix::mean() 
	{
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += element[i * columns + j];
			}
		}
		return answer / size();
	}

	real Matrix::sum() 
	{ 
		real answer = 0;  
		for (int b = 0; b < rows * columns; ++b) 
		{ 
			answer += element[b]; 
		} 
		return answer; 
	}

	real Matrix::sumSquares() 
	{ 
		real answer = 0;  
		for (int b = 0; b < rows * columns; ++b) 
		{ 
			answer += element[b] * element[b]; 
		} 
		return answer; 
	}

	real Matrix::geometricMean() {
		real answer = 1;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer *= element[i * columns + j];
			}
		}
		return pow(answer, real(1.0 / size()));
	}

	real Matrix::trace() 
	{
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += element[i * columns + i];
		}
		return answer;
	}

	int longestElement(Matrix& M)
	{
		int lelem = 0;
		for (int i = 0; i < M.rows; ++i)
		{
			for (int j = 0; j < M.columns; ++j)
			{
				int len = getLength(M(i,j));
				if (len > lelem)
					lelem = len;
			}
		}
		return lelem;
	}

	std::string Matrix::to_string(int precision) {
		std::ostringstream s;
		int defaultLength = 4;
		int largestElement = longestElement(*this);
		bool tf = isIntegerMatrix();
		if (tf) 
			--defaultLength;
		if (!tf) 
			++defaultLength;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				s <<
					std::fixed <<
					std::left <<
					std::setw(defaultLength + (precision) + (largestElement)) <<
					std::setprecision(precision) <<
					element[i * columns + j];
			}
			s << std::endl;
		}
		std::string temp = s.str();
		s.clear();
		return temp;
	}

	real Matrix::norm() 
	{//the Euclidean or Schur norm of the matrix, equal to the 
	 //square root of the sum of all squared elements.
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(element[i * columns + j], 2);
			}
		}
		return sqrt(answer);
	}

	real Matrix::FrobeniusNorm() 
	{ //the Froebenius norm is equal to the sqrt of the trace of 
	  //the matrix multiplied by its conjugate transpose.
		Matrix AT = *this;
		AT.transpose();
		return std::sqrt((*this * AT).trace());
	}

	real Matrix::pNorm(real p) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(element[i * columns + j], p);
			}
		}
		return sqrt(answer);
	}

	real Matrix::sumRow(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += element[r * columns + i];
		}
		return answer;
	}

	real Matrix::sumColumn(int r) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += element[i * columns + r];
		}
		return answer;
	}

	real Matrix::columnMean(int c) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += element[i * columns + c];
		}
		return (answer / rows);
	}

	real Matrix::rowMean(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += element[r * columns + i];
		}
		return (answer / columns);
	}

	real Matrix::rowNorm(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += pow(element[r * columns + i], 2);
		}
		return sqrt(answer);
	}

	real Matrix::columnNorm(int c) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[i * columns + c], 2);
		}
		return sqrt(answer);
	}

	real Matrix::columnSumSquares(int c) {
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[i * columns + c] - colMean, 2);
		}
		return answer;
	}

	real Matrix::rowSumSquares(int r) {
		real rMean = rowMean(r);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[r * columns + i] - rMean, 2);
		}
		return answer;
	}

	real Matrix::columnVariance(int c) {
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[i * columns + c] - colMean, 2);
		}
		return answer / (rows - 1);
	}

	real Matrix::columnCovariance(int c1, int c2) 
	{
		int n = rows;
		real answer = 0;
		Vector x = column(c1);
		Vector y = column(c2);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				answer += 0.5 * (x[i] - x[j]) * (y[i] - y[j]);
			}
		}
		return answer / (n * n);
	}

	real Matrix::columnDeviation(int c) 
	{
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[i * columns + c] - colMean, 2);
		}
		return sqrt(answer / (rows - 1));
	}

	real Matrix::rowCovariance(int r1, int r2) 
	{
		int n = columns;
		real answer = 0;
		Vector x = row(r1);
		Vector y = row(r2);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < n; ++j) {
				answer += 0.5 * (x[i] - x[j]) * (y[i] - y[j]);
			}
		}
		return answer / (n * n);
	}

	real Matrix::rowVariance(int r) {
		real rMean = rowMean(r);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[r * columns + i] - rMean, 2);
		}
		return answer / (columns - 1);
	}

	real Matrix::rowDeviation(int r) {
		real rMean = rowMean(r);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(element[r * columns + i] - rMean, 2);
		}
		return sqrt(answer / (columns - 1));
	}

	Vector Matrix::normalizedRow(int r) {
		Vector rw = row(r);
		rw *= rowNorm(r);
		return rw;
	}

	Vector Matrix::normalizedColumn(int c) {
		Vector col = row(c);
		col *= columnNorm(c);
		return col;
	}

	void Matrix::normalizeColumn(int c) {
		real temp = columnNorm(c);
		temp = 1 / temp;
		for (int i = 0; i < rows; ++i) {
			element[i * columns + c] *= temp;
		}
	}

	void Matrix::normalizeRow(int r) {
		real temp = rowNorm(r);
		temp = 1 / temp;
		for (int i = 0; i < columns; ++i) {
			element[r * columns + i] *= temp;
		}
	}

	void Matrix::normalize() {
		real d = norm();
		if (d != 0) {
			d = (1 / d);
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < columns; j++) {
					element[i * columns + j] *= d;
				}
			}
		}
	}

	real Matrix::sampleStandardDeviation() {
		real answer = 0;
		real Mean = mean();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(element[i * columns + j] - Mean, 2);
			}
		}
		return sqrt(answer / ((rows * columns) - 1));
	}

	real Matrix::populationStandardDeviation() {
		real answer = 0;
		real Mean = mean();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(element[i * columns + j] - Mean, 2);
			}
		}
		return sqrt(answer / (rows * columns));
	}

	real Matrix::ChiSquareTestStatistic() {
		real answer = 0;
		real df = (rows - 1) * (columns - 1);

		real sm = sum();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real Ei = (sumRow(i) * sumColumn(j) / sm);
				answer += pow(element[i * columns + j] - Ei, 2) / Ei;
			}
		}
		return answer;
	}

	real Matrix::ChiSquareDegreesFreedom() 
	{ 
		return  (rows - 1) * (columns - 1); 
	}

	void Matrix::removeRow(int n) {
		Vector newElements;
		newElements.resize((rows - 1) * columns, 0);
		int counter = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i != n) 
				{ 
					newElements[counter] = element[i * columns + j]; 
					++counter; 
				}
			}
		}
		element = newElements;
		--rows;
	}

	void Matrix::removeColumn(int n) {
		Vector newElements;
		newElements.resize(rows * (columns - 1), 0);
		int counter = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (j != n) 
				{ 
					newElements[counter] = element[i * columns + j]; 
					++counter; 
				}
			}
		}
		element = newElements;
		--columns;
	}

	void Matrix::trim() {//removes row and column vectors of 0's from inner, non-zero, matrix
		int firstNonZeroRow = 0;
		int firstNonZeroColumn = 0;

		//remove inner column/row 0 vectors
		for (int i = rows - 1; i >= 0; --i) {
			if (sumRow(i) == 0) { removeRow(i); }
			if (sumRow(i) != 0) { firstNonZeroRow = i; }
		}
		for (int i = columns - 1; i >= 0; --i) {
			if (sumColumn(i) == 0) { removeColumn(i); }
			if (sumColumn(i) != 0) { firstNonZeroColumn = i; }
		}

		//remove outer column/row 0 vectors;
		for (int i = 0; i < firstNonZeroRow; ++i) {
			if (sumRow(i) == 0) { removeRow(i); }
		}
		for (int i = 0; i < firstNonZeroColumn; ++i) {
			if (sumColumn(i) == 0) { removeColumn(i); }
		}
	}

	Matrix Matrix::directSum(Matrix A, Matrix B) {
		Matrix T = eye(A.rows + B.rows, A.columns + B.columns);
		for (int i = 0; i < T.rows; ++i) {
			for (int j = 0; j < T.columns; ++j) {
				if (i < A.rows && j < A.columns) 
					T(i, j) = A(i, j); 
				if (i >= A.rows && j >= A.columns) 
					T(i, j) = B(i - A.rows, j - A.columns);
			}
		}
		return T;
	}

	Matrix Matrix::tensorProduct(Matrix A, Matrix B) 
	{
		Matrix T = eye(A.rows * B.rows, A.columns * B.columns);
		int Rbound = (A.rows - 1);
		int Cbound = (A.columns - 1);
		int Rcounter = 0;
		int Ccounter = 0;
		for (int i = 0; i < T.rows; ++i) {
			for (int j = 0; j < T.columns; ++j) 
			{
				T(i, j) = A(Rcounter, Ccounter) * B(i % B.rows, j % B.columns);
				if (j % A.columns == Cbound && j) 
					++Ccounter;
			}
			Ccounter = 0;
			if (i % A.rows == Rbound && i) 
				++Rcounter;
		}
		return T;
	}

	bool canSwapOutZeroDiagonals(Matrix& M) {
		for (int i = 0; i < M.columns; ++i) 
		{
			if (!M.sumColumn(i)) 
				return false;
		}
		return true;
	}

	Matrix swapOutZeroDiagonals(Matrix M) 
	{
		Matrix M2 = M;
		while (zeroInDiag(M2)) 
		{
			for (int i = 0; i < M2.rows; ++i) {
				if (!M2(i, i)) 
				{
					for (int j = 0; j < M2.columns; ++j) 
					{
						if (M2(i, j) && M2(j, i) && i != j) 
						{   // Choose rows to swap to create M non-zero diagonal
							M2.swapRows(i, j); 
							j = M.columns;
						}
					}
				}
			}
		}
		return M2;
	}

	real Matrix::columnAbsMax(int n) {
		real m = element[n];
		if (rows <= 1)
			return m;
		for (int i = 0; i < rows; ++i) {
			real temp = element[i * columns + n];
			if (abs(temp) > m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnMax(int n) {
		real m = element[n];
		if (rows <= 1)
			return m;
		for (int i = 0; i < rows; ++i) {
			real temp = element[i * columns + n];
			if (temp > m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnAbsMin(int n) {
		real m = element[n];
		if (rows <= 1) 
			return m;
		for (int i = 1; i < rows; ++i) {
			real temp = element[i * columns + n];
			if (abs(temp) < m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnMin(int n) {
		real m = element[n];
		if (rows <= 1) { return m; }
		for (int i = 1; i < rows; ++i) {
			real temp = element[i * columns + n];
			if (temp < m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowAbsMax(int n) {
		real m = element[n * columns];
		if (columns <= 1) 
			return m;
		for (int i = 0; i < columns; ++i) {
			real temp = element[n * columns + i];
			if (abs(temp) > m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowMax(int n) {
		real m = element[n * columns];
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = element[n * columns + i];
			if (temp > m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowAbsMin(int n) {
		real m = element[n * columns];
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = element[n * columns + i];
			if (abs(temp) < m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowMin(int n) {
		real m = element[n * columns];
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = element[n * columns + i];
			if (temp < m) { m = temp; }
		}
		return m;
	}

	int Matrix::rowNonzeroValues(int rw) {
		int n = 0;
		for (int i = 0; i < columns; ++i) {
			if (std::abs(element[rw * columns + i]) > 0) {
				++n;
			}
		}
		return n;
	}

	int Matrix::columnNonzeroValues(int cl) {
		int n = 0;
		for (int i = 0; i < rows; ++i) {
			if (std::abs(element[i * columns + cl]) > 0) {
				++n;
			}
		}
		return n;
	}

	int Matrix::getPivot(int rw) {
		int val = 0;
		while (val < columns) {
			if (abs(element[rw * columns + val]) > 0) {
				return val;
			}
			++val;
		}
		return -1;
	}

	int Matrix::getReversePivot(int rw) {
		int val = columns - 1;
		while (val >= 0) {
			if (abs(element[rw * columns + val]) > 0) {
				return val;
			}
			--val;
		}
		return -1;
	}

	Matrix Matrix::GaussianElimination() 
	{
		int i, j, k;
		real piv;
		Matrix m(rows, columns, element);
		//Gauss-Jordan elimination
		for (i = 0; i < rows; i++) 
		{
			if (zeroInDiag(*this) && canSwapOutZeroDiagonals(*this)) 
			 //remove zeroes in diagonal if possible
				m = swapOutZeroDiagonals(m); 
			int ind = m.getPivot(i);
			if (ind >= 0) {
				piv = m(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) 
				{
					for (int l = i; l < rows; ++l) 
					{
						m(i, l) = m(i, l) / piv;
					}
					m(i, i) = 1;
				}
				//proceed down the column to make each non-pivot value = 0
				for (j = 0; j < rows; j++) 
				{
					if (m(j, i)) {
						if (j != i) {
							piv = m(j, i);
							for (k = i + 1; k < columns; k++) {
								m(j, k) = m(j, k) - (m(i, k) * piv);
							}
							m(j, i) = 0;
						}
					}
				}
			}
		}
		return m;
	}

	Matrix Matrix::lowerTriangularize() {
		int i, j, k;
		real piv;
		Matrix m(rows, columns, element);
		//Gauss-Jordan elimination
		for (i = rows - 1; i >= 0; --i) {
			if (zeroInDiag(*this) && canSwapOutZeroDiagonals(*this)) 
				m = swapOutZeroDiagonals(m);//remove zeroes in diagonal if possible
			int ind = m.getReversePivot(i);
			if (ind >= 0) {
				piv = m(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) {
					for (int l = 0; l < i; ++l) {
						m(i, l) = m(i, l) / piv;
					}
					m(i, i) = 1;
				}
				//proceed up the column to make each non-pivot value = 0
				for (j = i - 1; j >= 0; --j) {
					if (m(j, i) != 0) {
						if (j != i) {
							piv = m(j, i);//m(i,i);
							for (k = 0; k < i; ++k) {
								m(j, k) = m(j, k) - (m(i, k) * piv);
							}
							m(j, i) = 0;
						}
					}
				}
			}
		}
		return m;
	}

	Matrix Matrix::upperTriangularize() {
		int i, j, k;
		real piv;
		Matrix m(rows, columns, element);
		//Gauss-Jordan elimination
		for (i = 0; i < rows; i++) {
			if (zeroInDiag(*this) && canSwapOutZeroDiagonals(*this)) 
				m = swapOutZeroDiagonals(m); //remove zeroes in diagonal if possible
			int ind = m.getPivot(i);
			if (ind >= 0) {
				piv = m(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) {
					for (int l = i + 1; l < rows; ++l) {
						m(i, l) = m(i, l) / piv;
					}
					m(i, i) = 1;
				}
				//proceed down the column to make each non-pivot value = 0
				for (j = i + 1; j < rows; j++) {
					if (m(j, i) != 0) {
						if (j != i) {
							piv = m(j, i);//m(i,i);
							for (k = i + 1; k < columns; k++) {
								m(j, k) = m(j, k) - (m(i, k) * piv);
							}
							m(j, i) = 0;
						}
					}
				}
			}
		}
		return m;
	}

	int Matrix::rank() {
		int rowval = 0;
		for (int i = 0; i < rows; ++i) {
			real elem = sumRow(i);
			if (abs(elem) > 0) {
				++rowval;
			}
		}
		return rowval;
	}

	int Matrix::nullity() 
	{ 
		return columns - rank(); 
	}

	void Matrix::removeZeroColumns() 
	{
		for (int i = 0; i < columns; ++i) {
			if (!sumColumn(i))
			{
				removeColumn(i);
				--i;
			}
		}
	}

	void Matrix::removeZeroRows() 
	{
		for (int i = 0; i < rows; ++i) {
			if (!sumRow(i)) 
			{
				removeRow(i); 
				--i;
			}
		}
	}

	Vector Matrix::solve(Vector b)
	{//solve for Ax=b, where x will be the returned vector of variable values
	 //if m==n
		int n = rows;
		Vector x(n);
		for (int i = n - 1; i >= 0; i--) {
			x[i] = element[i * columns + n-1] / 
				element[i * columns + i];
			for (int k = i - 1; k >= 0; k--) {
				element[k * columns + n - 1] -= 
					element[k * columns + i] * x[i];
			}
		}
		return x;

	}

	Matrix Matrix::nullSpace() 
	{
		Matrix RRE = GaussianElimination();
		int kernelDimension = RRE.rank();//the rank of the matrix is equivalent to the dimension of the kernel
		//if (columns == kernelDimension) { return Matrix(1, 1); }

		//for a RRE matrix, the non-pivot columns are the free variables
		std::vector<int> piv = pivotColumns(RRE);
		RRE.removeZeroRows();

		std::vector<Vector> ans;
		int counter = 0;
		while (ans.size() < kernelDimension && counter < columns * 10) {
			for (int i = 0; i < piv.size(); ++i) {
				Matrix Mcopy(RRE);//set up matrix

				//make one of the pivot columns = 1, the others = 0. rotate which one is 1 over each iteration to find all distinct solution vectors
				for (int p = 0; p < piv.size(); ++p) 				
					if (p != i) 
						setColumnNonZeroValues(Mcopy, piv[p], 0);				
				setColumnNonZeroValues(Mcopy, piv[i], 1);

				//run the solver with for Ax=b where b is the 0 vector
				ans.push_back(Mcopy.solve(Vector(columns)));
				++counter;
			}
		}
		Matrix ansMat(ans);
		return ansMat;
	}

	Matrix Matrix::upperHessenbergForm() {
		if (rows != columns) { return Matrix(); }
		Matrix M(rows, columns, element);

		//elimination with 3 cases:  
		for (int j = 0; j < columns - 1; ++j) {
			real pivot = M(j + 1, j);
			for (int i = j + 2; i < rows; ++i) {
				Vector temp = M.row(j + 1);
				if (M(i, j) != 0) {
					if (abs(M(i, j)) == abs(pivot)) 
					{
						if (M(i, j) == pivot) 
						{ 
							temp*=-1; 
							goto skipUH; 
						}
						if ((M(i, j) == (-1 * pivot)) || ((-1 * M(i, j)) == pivot)) 
							goto skipUH;
					}
					if (abs(M(i, j)) != abs(pivot)) 					
						temp*=(-M(i, j) / pivot); 					

				skipUH:
					M.addToRow(i, temp);
				}
			}
		}
		return Matrix(rows, columns, M.element);
	}

	Matrix Matrix::lowerHessenbergForm() 
	{
		if (rows != columns)
			return Matrix();
		Matrix A(rows, columns, element);

		for (int j = columns - 1; j >= 0; --j) 
		{// Go the other way to remove the upper triangle of non-diagonal values.
			real pivot = A(j, j);
			for (int i = j - 1; i >= 0; --i) {
				if (i != j && A(i, j) != 0) {
					Vector temp = A.row(j);

					if (abs(A(i, j)) == abs(pivot)) {
						if (A(i, j) == pivot) {
							temp*= -1;
							goto skipLH;
						}
						if ((A(i, j) == (-1 * pivot)) || ((-1 * A(i, j)) == pivot)) 
						{ 
							goto skipLH; 
						}
					}
					if (abs(A(i, j)) != abs(pivot)) {
						real a = -A(i, j) / pivot;
						temp*= a;
					}

				skipLH:
					A.addToRow(i, temp);
				}
			}
		}

		return Matrix(rows, columns, A.element);
	}

	Matrix Matrix::characteristicMatrix() 
	{
		Matrix next(rows, columns, element);
		for (int k = 1; k < rows; ++k) {
			Matrix M = eye(columns, rows);
			for (int p = 1; p < columns + 1; ++p) 
			{
				if (p == (rows - k)) 
				{
					M.element[((rows - k - 1) * rows) + (p - 1)] = 
						(1 / (next.element[((rows - k) * columns) + (rows - k - 1)]));
				}
				else 
				{ 
					M.element[((rows - k - 1) * rows) + (p - 1)] = 
						-next.element[((rows - k) * columns) + (p - 1)] / 
						(next.element[((columns - k) * columns) + (columns - k - 1)]); 
				}
			}
			Matrix Minv = eye(columns, rows);
			for (int p = 1; p < columns + 1; ++p) 
			{ 
				Minv.element[(columns - k - 1) * columns + (p - 1)] = 
					next.element[(columns - k) * columns + (p - 1)]; 
			}
			next *= M;
			next = Minv * next;
		}

		for (int i = 1; i < next.rows; ++i) 
		{//blank out anything that's not in the first row;
			for (int j = 0; j < next.columns; ++j) {
				if (i == j) { next.element[i * columns + j] = 1; }
				else { next.element[i * columns + j] = 0; }
			}
		}
		return Matrix(rows, columns, next.element);
	}

	real Matrix::det() 	 
	{// Use Gauss-Jordan Elimination to get upper-triangle of matrix, then 
	 // multiply trace values to get the determinant.
		Matrix m = upperTriangularize();
		real det = 1;
		for (int i = 0; i < m.rows; ++i) 
			det *= m(i, i);		
		return det;
	}

	Matrix Matrix::reducedRowEchelonForm() 
	{
		Matrix A(rows, columns, element);

		while (zeroInDiag(A)) 
		{
			for (int i = 0; i < rows; ++i) 
			{
				if (A(i, i) == 0) 
				{
					for (int j = 0; j < columns; ++j) 
					{// Choose rows to swap to create a non-zero diagonal
						if (A(i, j) != 0 && A(j, i) != 0 && i != j) 
						{ 
							A.swapRows(i, j);
							j = columns;
						}
					}
				}
			}
		}

		// Iterate and eliminate the lower triangle of values
		for (int j = 0; j < rows; ++j) 
		{
			real pivot = A(j, j);
			for (int i = 0; i < rows; ++i) 
			{
				if (i != j && A(i, j) != 0) 
				{
					Vector temp = A.row(j);
					if (abs(A(i, j)) == abs(pivot)) 
					{
						if (A(i, j) == pivot) 
						{
							temp*= -1;
							goto skip;
						}
						if ((A(i, j) == (-1 * pivot)) || ((-1 * A(i, j)) == pivot)) 
							goto skip;
					}
					if (abs(A(i, j)) != abs(pivot)) 
					{
						real a = -A(i, j) / pivot;
						temp*= a;
					}

				skip:
					A.addToRow(i, temp);

				}
			}
		}

		for (int i = 0; i < rows; ++i) 
		{//make diagonals into all 1's
			if (A(i, i) != 0 && A(i, i) != 1) {
				real a = ((real)1 / (A(i, i)));
				A.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, A.element);
	}

	Matrix Matrix::addRow(Vector vec) 
	{//adds vector to right hand side by default
		Matrix M(rows + 1, columns);
		int rw = rows;
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (i == rw) {
					M(i, j) = vec[j];
				}
				else {
					if (i < rw) {
						M(i, j) = element[i * columns + j];
					}
					if (i > rw) {
						M(i, j) = element[(i - 1) * columns + j];
					}
				}
			}
		}
		return M;
	}

	Matrix Matrix::addRow(int rw, Vector vec) 
	{
		Matrix M(rows + 1, columns);
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (i == rw) {
					M(i, j) = vec[j];
				}
				else {
					if (i < rw)
						M(i, j) = element[i * columns + j];					
					if (i > rw) 
						M(i, j) = element[(i-1) * columns + j];					
				}
			}
		}
		return M;
	}

	Matrix Matrix::addColumn(Vector vec) 
	{//adds to bottom of matrix by default.
		Matrix M(rows, columns + 1);
		int cl = columns;
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (j == cl) {
					M(i, j) = vec[i];
				}
				else {
					if (j < cl)
						M(i, j) = element[i * columns + j];					
					if (j > cl)
						M(i, j) = element[i * columns + (j-1)];					
				}
			}
		}
		return M;
	}

	Matrix Matrix::addColumn(int cl, Vector vec) 
	{
		Matrix M(rows, columns + 1);
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (j == cl) {
					M(i, j) = vec[i];
				}
				else {
					if (j < cl)
						M(i, j) = element[i * columns + j];
					if (j > cl) 
						M(i, j) = element[i * columns + (j - 1)];
				}
			}
		}
		return M;
	}

	bool Matrix::isInconsistent(Vector b) 
	{//by Rouche-Capelli Theorem, any system of linear equations 
	 //(underdetermined or otherwise) is inconsistent if the rank of the 
	 //augmented matrix is greater than the rank of the coefficient matrix.
		if (rank() <= addColumn(columns - 1, b).rank()) 
			return true;
		return false;
	}

	std::vector<Matrix> Matrix::SVD() {
		/*Takes an mxn matrix a and decomposes it into udv, where u,v are
		left and right orthogonal transformation matrices, and d is a
		diagonal matrix of singular values.*/
		std::vector<Matrix> Mat;
		Vector w(columns);
		Vector eigens1, eigens2;
		Matrix v(rows, columns);
		Matrix M(rows, columns, element);
		Matrix MT = M.transpose();

		//search for U
		Matrix MMT = M * MT;
		Mat.push_back(MMT.eigenvectors());

		//search for V
		Matrix MTM = MT * M;
		Mat.push_back(MTM.eigenvectors());

		//search for S
		eigens1 = MMT.eigenvaluesRealExact();
		eigens2 = MTM.eigenvaluesRealExact();
		eigens1.insert(eigens1.begin(), eigens2.begin(), eigens2.end());
		eigens2.clear();
		eigens1.removeNullValues();
		Mat.push_back(Matrix(eigens1));
		return(Mat);
	}

	std::vector<Matrix> LU(Matrix& M) 
	{
		std::vector<Matrix> answer;
		Matrix A = M;
		Matrix L(M.rows, M.columns);
		Matrix U(M.rows, M.columns);

		int i, j, k;
		real sum = 0;

		// Set row diagonals to 1.
		for (i = 0; i < M.rows; i++)
			U(i, i) = 1;

		for (j = 0; j < M.columns; j++) {
			for (i = j; i < M.rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L(i, k) * U(k, j);
				}
				L(i, j) = A(i, j) - sum;
			}

			for (i = j; i < M.rows; i++) 
			{
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L(j, k) * U(k, i);
				}
				real divisor = L(j, j);
				if (!divisor)//prevent division by 0.
					divisor = std::numeric_limits<real>::min();
				U(j, i) = (A(j, i) - sum) / divisor;
			}
		}
		answer.push_back(L);
		answer.push_back(U);
		return answer;
	}

	std::vector<Matrix> Matrix::QR() 
	{//QR decomposition by Gram-Schmidt process
		std::vector<Matrix> A;
		Matrix Q(rows, columns);
		Matrix R(rows, columns);

		//set up Q matrix
		Q.setColumn(0, column(0).length());
		for (int i = 0; i < columns; ++i) {
			Vector c = column(i);
			Vector a = c;
			for (int j = i; j > 0; --j) {
				Vector e = Q.column(j - 1);
				real dp = dot(c, e);
				Vector tp = e*dp;
				a -= tp;
			}
			a = length(a);
			Q.setColumn(i, a);
		}

		//set up R matrix
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (j >= i) {
					Vector a = column(j);
					Vector b = Q.column(i);
					real temp = dot(a, b);
					R(i, j) = temp;
				}
			}
		}
		A.push_back(Q);
		A.push_back(R);
		return A;
	}

	Matrix Matrix::inverseByQR() 
	{
		if (!det()) 
			return Matrix();
		std::vector<Matrix> qr = QR();
		Matrix qr_t = qr[0].transpose();
		return qr[1] * qr_t;
	}

	Matrix Matrix::GramSchmidt() 
	{	//Orthogonalization by Gram-Schmidt process using Gaussian Elimination.
		if (!det())
			return Matrix();
		Matrix Aug(rows, columns, element);
		Matrix AT(rows, columns, element);
		AT = AT.transpose();
		Matrix A = *this * AT;

		while (zeroInDiag(A)) 
		{
			for (int i = 0; i < rows; ++i) 
			{
				if (A(i, i)) 
				{
					for (int j = 0; j < columns; ++j) 
					{
						//choose rows to swap to create a non-zero diagonal
						if (A(i, j) != 0 && A(j, i) != 0 && i != j) 
						{
							A.swapRows(i, j);
							Aug.swapRows(i, j);
							j = columns;
						}
					}
				}
			}
		}

		for (int j = 0; j < columns; ++j) 
		{	// Iterate down matrix and eliminate the lower triangle of values.
			real pivot = A(j, j);
			for (int i = j + 1; i < rows; ++i) {
				if (i != j && A(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A(i, j)) == abs(pivot)) {
						if (A(i, j) == pivot) {
							temp *= -1;
							temp2 *= -1;
							goto skip;
						}
						if ((A(i, j) == (-1 * pivot)) || ((-1 * A(i, j)) == pivot)) 
							goto skip;
					}
					if (abs(A(i, j)) != abs(pivot)) 
					{
						real a = -A(i, j) / pivot;
						temp *= a;
						temp2 *= a;
					}

				skip:
					A.addToRow(i, temp);
					Aug.addToRow(i, temp2);

				}
			}
		}

		for (int i = 0; i < rows; ++i) 
		{	//Make diagonals into all 1's.
			if (A(i, i) != 0 && A(i, i) != 1) {
				real a = ((real)1 / (A(i, i)));
				A.multiplyRow(i, a);
				Aug.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, Aug.element);
	}

	std::vector<int> JacobiIndexing(Matrix& M) 
	{// Returns indices i,j for largest element in matrix quickly for use in 
	 //'JacobiTransformation()' function.
		std::vector<int> answer;
		answer.push_back(-1);
		answer.push_back(-1);
		real largest = 0;
		for (int i = 0; i < M.rows; ++i) 
		{
			for (int j = 0; j < M.columns; ++j) 
			{
				if (i != j && std::abs(M(i,j)) > largest) 
				{
					largest = std::abs(M(i, j));
					answer[0] = i;
					answer[1] = j;
				}
			}
		}
		return answer;
	}

	Matrix Matrix::GivensRotationMatrix(real a1, real a2, int r1, int r2) 
	{//eliminate a2 with a1, a1 is in row r1, a2 is in row r2.
		Matrix el = eye(rows, rows);

		real r = sqrt((a1 * a1) + (a2 * a2));
		real s = a1 / r;
		real negs = -s;
		real cs = a2 / r;

		el.element[r1 * rows + r1] = cs;
		el.element[r2 * rows + r2] = cs;
		el.element[r1 * rows + r2] = s;
		el.element[r2 * rows + r1] = negs;
		return Matrix(rows, rows, el.element);
	}

	Matrix JacobiRotationMatrix(Matrix& M, int p, int q, real c, real s) 
	{	//Helper function used in Jacobi Transformation.
		Matrix el = eye(M.rows, M.rows);
		el.element[p * M.rows + p] = c;
		el.element[q * M.rows + q] = c;
		el.element[p * M.rows + q] = s;
		el.element[q * M.rows + p] = -s;
		return Matrix(M.rows, M.rows, el.element);
	}

	std::vector<Matrix> JacobiTransformation(Matrix& M, highpUint limit)
	{	// Outputs an array with both diagonalized matrix and eigenvectors
		Matrix D = M;//matrix to rotate/diagonalize
		Matrix Eigenvectors = eye(M.rows, M.columns);//matrix of Eigenvectors produced by successive Jacobi rotations.

		int counter = 0;//counts the number of iterations, to stop at some set limit.	
		while (!D.isDiagonalized() && counter < limit) 
		{
			std::vector<int> elem = JacobiIndexing(D);
			if (elem[0] > D.rows || elem[1] > D.columns || elem[0] < 0 || elem[1] 
				< 0 || elem[0] == elem[1]) 
				break;
			int p = elem[0];
			int q = elem[1];
			real angle = (D(q, q) - D(p, p)) / (2 * D(p, q));
			int sign = 1;
			if (angle < 0) 
				sign = -1;
			real t = sign / (abs(angle) + sqrt((angle * angle) + 1));
			real c = 1 / sqrt((t * t) + 1);
			real s = t * c;
			Matrix rot = JacobiRotationMatrix(D, p, q, c, s);
			Matrix rotT = rot.transpose();
			D = rotT * D;
			D *= rot;
			Eigenvectors = Eigenvectors * rot;
			++counter;
		}
		std::vector<Matrix> answer;
		answer.push_back(D);
		answer.push_back(Eigenvectors);
		return answer;
	}

	Matrix Matrix::eigenvectors(int iterations)
	{
		if (isSymmetric()) //If a symmetric matrix, simply use Jacobi transformation method.
		{
			std::vector<Matrix> EV = JacobiTransformation(*this);
			return EV[1];
		}
		else //Else, use iterative QR decomp method.
		{
			std::vector<Matrix> qr = QR();
			Matrix Eigenvalues = qr[0];
			Matrix A2 = qr[1] * qr[0];

			for (int i = 1; i < iterations; ++i) 
			{	// Then, QR decomposition. Use Q as multiplicand where 
				// Eigenvectors = Q1*Q2*...Qn
				qr = A2.QR();
				if (qr[0].rows == 0 && qr[1].rows == 0) 
					break;
				Eigenvalues = Eigenvalues * qr[0];
				A2 = qr[1] * qr[0];
				qr.clear();
			}
			return Eigenvalues;
		}
		return Matrix();
	}

	Matrix Matrix::Cholesky() 
	{	//returns matrix L decomposition of matrix where A = LL^T.
		if (rows <= 1 || columns <= 1) 
			return Matrix();
		Matrix L(rows, columns);
		L(0, 0) = std::sqrt(element[0]);

		for (int j = 0; j < columns; ++j) 
		{
			for (int i = j; i < rows; ++i) 
			{
				if (i == j) 
				{
					real temp = element[j * columns + j];
					real temp2 = 0;
					for (int k = 0; k < j; ++k) {
						temp2 += pow(L(j, k), 2);
					}
					L(j, j) = std::sqrt(temp - temp2);
				}

				if (i > j) {
					real temp = (1 / L(j, j));
					real temp2 = element[i * columns + j];
					real temp3 = 0;
					for (int k = 0; k < j; ++k)
						temp3 += (L(i, k) * L(j, k));

					if (temp) 
						L(i, j) = temp * (temp2 - temp3);
				}
			}
		}
		return Matrix(rows, columns, L.element);
	}

	std::vector<Matrix> Matrix::LDL() 
	{	//returns matrix L and D decomposition of matrix where A = LDL^T.
		std::vector<Matrix> answer;
		if (rows <= 1 || columns <= 1)
			return answer;
		Matrix L(rows, columns);
		Matrix D(rows, columns);
		L(0, 0) = std::sqrt(element[0]);

		for (int j = 0; j < columns; ++j) {
			for (int i = j; i < rows; ++i) {
				if (i == j) {
					L(j, j) = 1;
					real temp = 0;
					for (int k = 0; k < j; ++k) {
						temp += pow(L(j, k), 2) * D(k, k);
					}
					D(j, j) = element[j * columns + j] - temp;
				}

				if (i > j) {
					real temp = (1 / D(j, j));
					real temp2 = element[i * columns + j];
					real temp3 = 0;
					for (int k = 0; k < j; ++k) {
						temp3 += (L(i, k) * L(j, k) * D(k, k));
					}
					L(i, j) = temp * (temp2 - temp3);
				}
			}
		}

		answer.push_back(L);
		answer.push_back(D);
		return answer;
	}

	Matrix Matrix::dominantEigenvector(int iterations) 
	{//calculated by power method
		Matrix A(rows, columns, element);
		Matrix rd(A.rows, 1);
		rd(0, 0) = 1;//set first value of matrix to 1
		Matrix init = rd;
		for (int i = 0; i < iterations; ++i) {
			rd = A * rd;
			rd = rd * (1.0 / rd.element[A.rows - 1]);
		}
		return rd;//dominant eigenvector
	}

	real Matrix::dominantEigenvalue(int iterations) 
	{//calculated by power method
		Matrix A = *this;
		Matrix rd(A.rows, 1);
		rd(0, 0) = 1;//set first value of matrix to 1
		Matrix init = rd;
		for (int i = 0; i < iterations; ++i)
		{
			rd = A * rd;
			rd = rd * (1.0 / rd.element[A.rows - 1]);
		}
		Matrix step1 = A * rd;
		step1 *= rd;
		real egv = step1.sumColumn(0);
		egv /= (rd * rd).sumColumn(0);
		return egv;//dominant eigenvalue
	}

	bool Matrix::isPositiveDefinite() 
	{//By theorem, if all eigenvalues are positive, the matrix is positive definite.
		Vector vec = eigenvaluesRealExact();
		for (int i = 0; i < vec.size(); ++i)
			if (vec[i] <= 0) 
				return false;		
		return true;
	}

	bool Matrix::isPositiveSemidefinite() 
	{//if all eigenvalues are positive or 0, the matrix is positive semi-definite
		Vector vec = eigenvaluesRealExact();
		for (int i = 0; i < vec.size(); ++i) {
			if (vec[i] < 0) { return false; }
		}
		return true;
	}

	Vector Matrix::eigenvaluesByGaussianElimination() {
		addColumn(Vector(rows));//make augmented matrix w/ column vector of 0's for solution
		Matrix GE = GaussianElimination();
		return GE.characteristicPolynomial().realRoots();
	}

	Vector Matrix::eigenvaluesRealExact() 
	{//finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step
		Polynomial p = characteristicPolynomial();
		return p.realRoots();
	}

	std::vector<ComplexNumber> Matrix::eigenvaluesRealAndComplexExact(int iterations) 
	{	//Finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step,
		//then either factors those roots to get the last few complex roots or uses the 
		//numerical Durand-Kerner method of root-finding.

		Polynomial p = characteristicPolynomial();
		Vector realRoots = p.realRoots();
		std::vector<ComplexNumber> complexRoots;
		for (int i = 0; i < realRoots.size(); ++i) {
			//std::wcout << realRoots[i] << L"\n";//debug
			complexRoots.push_back(ComplexNumber(realRoots[i], 0));
		}

		Vector bounds = p.getComplexRootBounds();
		real lowbound = bounds[0];
		real highbound = bounds[1];
		real lowboundC = bounds[0];

		for (int i = 0; i < realRoots.size(); ++i) {
			if (std::abs(int(p.size()) - (int)realRoots.size()) <= 3) 
			{// Use the shift method only if the real roots reduce it down to an easily 
			 // calculable size (order 2 or less), otherwise the results will be inaccurate.
				if (p.size() <= 3) 
				{ 
					i = realRoots.size(); 
					break; 
				}
				if (p.size() > 3) 				 
					p = p.factorOutBinomial(realRoots[i]); 				
			}
		}
		if (p.size() < 3) 
			return complexRoots;

		if (p.size() == 3) 
		{
			real part2 = pow(p[1], 2) - (4 * p[0] * p[2]);
			if (part2 < 0) {
				part2 = part2 / (2 * p[2]);
				real part1 = (-p[1] / (2 * p[2]));
				complexRoots.push_back(ComplexNumber(part1, part2));
				complexRoots.push_back(ComplexNumber(part1, -part2));
			}
			return complexRoots;
		}

		if (p.size() > 3) 
		{	// Start finding roots with Durand-Kerner method.
			ComplexNumber z = ComplexNumber(lowbound, lowboundC);
			int size = sizeof(z);
			std::vector<ComplexNumber> R;
			for (int i = 0; i < p.size(); i++)
				R.push_back(std::pow(z, i)); 			

			for (int i = 0; i < iterations; i++) {
				for (int j = 0; j < p.size(); j++) {
					ComplexNumber B = p.evaluate(R[j]);
					for (int k = 0; k < p.size(); k++) {
						if (k != j) { B /= (R[j] - R[k]); }
					}
					R[j] -= B;
				}
			}

			for (int i = 0; i < R.size(); ++i) 
			{	// Now filter out all the bad results.
				if (R[i].im() != 0) 
				{	// Only complex roots accepted.
					R[i] = p.NewtonsMethodComplex(p,R[i]);
					ComplexNumber temp = p.evaluate(R[i]);
					real a = std::abs(temp.re());
					real b = std::abs(temp.im());
					if (a < 0.1 && b < 0.1) {
						a = std::abs(R[i].re());
						b = std::abs(R[i].im());
						bool isSimilar = false;
						for (int i = 0; i < complexRoots.size(); ++i) {
							real x = std::abs(complexRoots[i].re());
							real y = std::abs(complexRoots[i].im());
							if (std::abs(a - x) < 0.001 && std::abs(b - y) < 0.001) 
								isSimilar = true;
						}
						if (!isSimilar) 
						{	// If this is indeed a new root, save the root and its conjugate.
							complexRoots.push_back(R[i]);
							complexRoots.push_back(ComplexNumber(R[i].re(), R[i].im() * -1));
						}
					}
				}
			}
			R.clear();
		}
		return complexRoots;
	}

	Vector Matrix::eigenvaluesNumerical() 
	{//calculated by power method w/shifts
		Matrix A = *this;
		Matrix v;
		Vector egv;
		for (int i = 0; i < A.rows; ++i) 
		{	//A' = A - (eigenvalue)*U*UT
			v = A.dominantEigenvector();
			Matrix A_1 = A * v;
			A_1 *= v;
			real eigenval = A_1.sumColumn(0);
			real val = (v * v).sumColumn(0);
			if (val != 0) { eigenval /= val; }
			else { eigenval = 0; }
			egv.push_back(eigenval);
			v.normalize();

			Matrix v_t = v.transpose();
			Matrix U = v * v_t;

			U = U * (-1.0 * eigenval);
			A = A.add(A, U);
		}
		return egv;
	}

	Matrix Matrix::eigenvalueMatrix() 
	{
		Matrix A = *this;
		std::vector<Matrix> QR;

		int iterations = pow((A.rows * A.columns), 2);
		for (int i = 0; i < iterations; ++i) 
		{//first, LU decomposition
			QR = LU(A);
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}

		for (int i = 0; i < iterations; ++i) 
		{//then, QR decomposition
			QR = A.QR();
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}
		for (int i = 0; i < iterations; ++i) 
		{//then, QR decomposition
			QR = A.HouseholderQR();
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}

		return A;
	}

	Matrix Matrix::eigenvalueMatrix(Matrix A, int loops) 
	{
		if (loops <= 1) 		
			return A; 		
		--loops;
		std::vector<Matrix> QR;

		int iterations = (A.rows * A.columns);
		for (int i = 0; i < iterations; ++i) 
		{//first, LU decomposition
			QR = LU(A);
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}

		for (int i = 0; i < iterations; ++i) 
		{//then, QR decomposition
			QR = A.QR();
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}
		for (int i = 0; i < iterations; ++i) 
		{//then, QR decomposition
			QR = A.HouseholderQR();
			if (QR[0].rows == 0 && QR[1].rows == 0) 
				break;
			A = QR[1] * QR[0];
		}

		return eigenvalueMatrix(A, loops);
	}

	Matrix Matrix::leastSquares(Matrix X, Matrix Y) 
	{//B=(X^T*X)^-1 * X^T * Y  
		Matrix XT = X.transpose();
		XT *= X;
		Matrix XTinv = XT.inverse();
		XTinv *= XT;
		return XTinv * Y;
	}

	Matrix Matrix::leastSquaresMatrix() 
	{//B=(X^T*X)^-1 * X^T * Y 
	 //NOTE: the matrix must be formatted in MINITAB format, so that the first column 
	 //contains the Y values, the rest are the column vectors of the X_i variable values.

		Matrix M(rows, columns, element);//copy parent matrix
		if (columns != 2 && rows == 2) 
			M = M.transpose(); //put matrix into column vector form

		Matrix X(M.rows, M.columns, element);//isolate X variable values
		for (int i = 0; i < X.rows; ++i) 
			X(i, 0) = 1;		

		Matrix Y(M.rows, 1, M.column(0));//create seperate vector of y values
		Matrix lhs = X.transpose();
		lhs *= X;
		lhs = lhs.inverse();
		Matrix rhs = X.transpose();
		rhs *= Y;
		return lhs * rhs;
	}

	Function Matrix::leastSquares() 
	{	//B=(X^T*X)^-1 * X^T * Y 
		//NOTE:  the matrix must be formatted in MINITAB format, so that the first column contains
		//the Y values, the rest are the column vectors of the X_i variable values. 
		//See https://onlinecourses.science.psu.edu/stat501/node/292 for examples. 

		Matrix coef = leastSquaresMatrix();
		Vector cf;
		std::string funct = "f(";
		for (int i = 0; i < coef.rows; ++i) {
			cf.push_back(coef(i, 0));
			funct += (i+'a') + ",";
		}
		funct.pop_back();
		funct += ") = ";

		for (int i = 0; i < cf.size(); ++i) {
			funct += to_string_precision(cf[i],30);
			if (i > 0) { funct += "*" + (i - 1 + 'a'); }
			funct += " + ";
		}
		//Remove final "+ " chars from string
		funct.pop_back();
		funct.pop_back();

		return Function(funct);
	}

	Vector Matrix::residuals() 
	{//get least squares regression function
		Function f = leastSquares();

		//get y-hat values
		Vector yhat;
		for (int i = 0; i < rows; ++i) {
			Vector rw = row(i);
			Vector vars;
			for (int j = 1; j < columns; ++j) {
				vars.push_back(rw[j]);
			}
			yhat.push_back(f.evaluate(vars));
		}

		// Calculate resid. = y-yhat. In MINITAB format, column 0 is the y-value column.
		Vector residuals;
		for (int i = 0; i < yhat.size(); ++i)
			residuals.push_back(element[i * columns] - yhat[i]);

		return residuals;
	}

	real Matrix::Rsquared() 
	{// aka Pearsson's "R"^2. Ref: https://people.richland.edu/james/ictcm/2004/multiple.html

		Function f = leastSquares();//get leastSquares multivariate function

		//get y-hat values
		Vector yhat;
		for (int i = 0; i < rows; ++i) {
			Vector rw = row(i);
			Vector vars;
			for (int j = 1; j < columns; ++j) {
				vars.push_back(rw[j]);
			}
			yhat.push_back(f.evaluate(vars));
		}

		//get residuals and residual standard deviation
		Vector resi = residuals();
		real residualSD = resi.sampleStandardDeviation();

		real ymean = columnMean(0);	//ymean

		//fill arrays of error calculations
		Vector explainedError;
		Vector unexplainedError;
		Vector totalError;
		for (int i = 0; i < resi.size(); ++i) 
		{
			explainedError.push_back(pow(yhat[i] - ymean, 2));// sum(yhat - yman)^2
			unexplainedError.push_back(pow(element[i * columns] - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(element[i * columns] - ymean, 2));// sum(y - ymean)^2
		}

		//degress of freedom calculations
		real dfRegression = columns - 1;	//#predictor varibles aka 'k'
		real dfTotal = rows - 1;//total degrees freedom
		real dfResiduals = dfTotal - dfRegression;//# of residuals

		//sum of squares calculations
		real SSregression = explainedError.sum();
		real SSresiduals = unexplainedError.sum();
		real SStotal = totalError.sum();

		//mean square (MS) calculations
		real MSregression = SSregression / dfRegression;
		real MSresiduals = SSresiduals / dfResiduals;
		real MStotal = SStotal / dfTotal;

		//important statistics
		real sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
		real FTestStatistic = MSregression / MSresiduals;
		real Rsquared = (SStotal - SSresiduals) / SStotal;
		real adjustedRsquared = (MStotal - MSresiduals) / MStotal;

		return Rsquared;
	}

	real Matrix::adjustedRsquared() 
	{// aka Pearsson's "R"^2. Ref: https://people.richland.edu/james/ictcm/2004/multiple.html

		Function f = leastSquares();//get leastSquares multivariate function

		//get y-hat values
		Vector yhat;
		for (int i = 0; i < rows; ++i) {
			Vector rw = row(i);
			Vector vars;
			for (int j = 1; j < columns; ++j) {
				vars.push_back(rw[j]);
			}
			yhat.push_back(f.evaluate(vars));
		}

		//get residuals and residual standard deviation
		Vector resi = residuals();
		real residualSD = resi.sampleStandardDeviation();

		real ymean = columnMean(0);	//ymean

		//fill arrays of error calculations
		Vector explainedError;
		Vector unexplainedError;
		Vector totalError;
		for (int i = 0; i < resi.size(); ++i) {
			explainedError.push_back(pow(yhat[i] - ymean, 2));// sum(yhat - yman)^2
			unexplainedError.push_back(pow(element[i * columns] - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(element[i * columns] - ymean, 2));	// sum(y - ymean)^2
		}

		//degress of freedom calculations
		real dfRegression = columns - 1;	//#predictor varibles aka 'k'
		real dfTotal = rows - 1;//total degrees freedom
		real dfResiduals = dfTotal - dfRegression;//# of residuals

		//sum of squares calculations
		real SSregression = explainedError.sum();
		real SSresiduals = unexplainedError.sum();
		real SStotal = totalError.sum();

		//mean square (MS) calculations
		real MSregression = SSregression / dfRegression;
		real MSresiduals = SSresiduals / dfResiduals;
		real MStotal = SStotal / dfTotal;

		//important statistics
		real sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
		real FTestStatistic = MSregression / MSresiduals;
		real Rsquared = (SStotal - SSresiduals) / SStotal;
		real adjustedRsquared = (MStotal - MSresiduals) / MStotal;

		return adjustedRsquared;
	}

	real Matrix::FTestStatistic() 
	{//Ref: https://people.richland.edu/james/ictcm/2004/multiple.html

		Function f = leastSquares();//get leastSquares multivariate function

		//get y-hat values
		Vector yhat;
		for (int i = 0; i < rows; ++i) {
			Vector rw = row(i);
			Vector vars;
			for (int j = 1; j < columns; ++j) {
				vars.push_back(rw[j]);
			}
			yhat.push_back(f.evaluate(vars));
		}

		//get residuals and residual standard deviation
		Vector resi = residuals();
		real residualSD = resi.sampleStandardDeviation();

		real ymean = columnMean(0);	//ymean

		//fill arrays of error calculations
		Vector explainedError;
		Vector unexplainedError;
		Vector totalError;
		for (int i = 0; i < resi.size(); ++i) {
			explainedError.push_back(pow(yhat[i] - ymean, 2));// sum(yhat - yman)^2
			unexplainedError.push_back(pow(element[i * columns] - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(element[i * columns] - ymean, 2));// sum(y - ymean)^2
		}

		//degress of freedom calculations
		real dfRegression = columns - 1;//#predictor varibles aka 'k'
		real dfTotal = rows - 1;//total degrees freedom
		real dfResiduals = dfTotal - dfRegression;//# of residuals

		//sum of squares calculations
		real SSregression = explainedError.sum();
		real SSresiduals = unexplainedError.sum();
		real SStotal = totalError.sum();

		//mean square (MS) calculations
		real MSregression = SSregression / dfRegression;
		real MSresiduals = SSresiduals / dfResiduals;
		real MStotal = SStotal / dfTotal;

		//important statistics
		real sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
		real FTestStatistic = MSregression / MSresiduals;
		real Rsquared = (SStotal - SSresiduals) / SStotal;
		real adjustedRsquared = (MStotal - MSresiduals) / MStotal;

		return FTestStatistic;
	}

	std::vector<Matrix> Matrix::thinQR() 
	{
		std::vector<Matrix> answer = QR();
		answer[0];
		answer[1].trim();
		return answer;
	}

	std::vector<Matrix> Matrix::HouseholderQR() 
	{//given some matrix, this produces its Householder QR factorization
		Matrix A = *this;//copy this matrix.
		Matrix Q = eye(rows, columns);
		for (int j = 0; j < columns - 1; ++j) {
			Vector temp = A.column(j);
			Vector vec;
			for (int i = j; i < rows; ++i) { vec[i - j] = temp[i]; }
			real temp1 = sgn(vec[0]) * vec.length();
			vec[0] += temp1;
			vec = vec / vec.length();//normalize.
			Matrix Vec(rows - j, 1, vec);
			Matrix H = Vec.GramMatrix(Vec);
			H *= 2;
			Matrix I = eye(rows - j, columns - j);
			H = I - H;
			if (j > 0) { H = H.expandToUpperLeft(rows, columns); }
			A = H * A;
			Q *= H;
		}
		std::vector<Matrix> QR;
		QR.push_back(Q);
		QR.push_back(A);
		return QR;
	}

	std::vector<Matrix> Matrix::thinHouseholderQR() {
		std::vector<Matrix> answer = HouseholderQR();
		answer[0];
		answer[1].trim();
		return answer;
	}

	Matrix Matrix::pseudoinverse()
	{//calculates the Moore-Penrose inverse: A* = R^-1 * Q^T.  Note:  A*A = I.
		std::vector<Matrix> qr = thinQR();
		Matrix inv1 = qr[1].inverse();
		Matrix inv2 = qr[0].transpose();
		Matrix result = inv1 * inv2;
		return result;
	}

	Polynomial Matrix::characteristicPolynomial() {
		Vector ans = characteristicMatrix().row(0);
		Vector answer(columns + 1);
		for (int i = 0; i < columns; ++i) 
			answer[columns - 1 - i] = -ans[i];
		answer[columns] = 1;
		return Polynomial(answer);
	}

	Matrix directionMatrix(Vector v) 
	{
		int sz = v.size() + 1;
		Matrix A = eye(sz, sz);
		for (int i = 0; i < v.size(); ++i) 
			A(i, i) = v[i];
		return A;
	}

	Matrix eye(int n)
	{ 
		return eye(n, n);
	}

	Matrix eye(int n, int m)
	{
		Vector vals(n * m);
		for (int i = 0; i < n; ++i) 
		{
			for (int j = 0; j < m; ++j) 
			{
				vals[i * m + j] = 0;
				if (i == j) 
					vals[i * m + j] = 1;
			}
		}
		return Matrix(n, m, vals);
	}

	Matrix positionMatrix(Vector v) 
	{
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		for (int i = 0; i < v.size(); ++i) 
			A(i, i) = v[i];
		return A;
	}

	Matrix rotationMatrix(Vector v) {
		//Input rotation vector should be of form <angle,x,y,z> for 3D rotation, 
		//or <angle,x,y> for 2D. If x, y, or z == 1, then it is a rotation around 
		//that axis.  Only angle and one coordinate should be non-zero.
		if (v.size() < 3) 
			return Matrix();
		if (v.size() == 3) 
		{//2D rotation matrix
			Matrix M(2, 2);
			M(0, 0) = (v[0]);
			M(0, 1) = (v[0]);
			M(1, 0) = (v[0]);
			M(1, 1) = (v[0]);
			return M;
		}
		if (v.size() == 4) 
		{//3D rotation matrix (contained in a 4x4 homogeneous coordinate matrix)
			if (v[1] == 1) 
			{//x-axis rotation
				Matrix M = eye(4, 4);
				M(1, 1) = (v[0]);
				M(1, 2) = (v[0]);
				M(2, 1) = (v[0]);
				M(2, 2) = (v[0]);
				return M;
			}
			if (v[2] == 1) 
			{//y-axis rotation
				Matrix M = eye(4, 4);
				M(0, 0) = (v[0]);
				M(0, 2) = (v[0]);
				M(2, 0) = (v[0]);
				M(2, 2) = (v[0]);
				return M;
			}
			if (v[3] == 1) 
			{//z-axis rotation
				Matrix M = eye(4, 4);
				M(0, 0) = (v[0]);
				M(0, 1) = (v[0]);
				M(1, 0) = (v[0]);
				M(1, 1) = (v[0]);
				return M;
			}
		}
		return Matrix();
	}

	Matrix scalingMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A = eye(sz, sz);
		for (int i = 0; i < v.size(); ++i) 
			A(i, i) = v[i];
		return A;
	}

	Matrix translationMatrix(Vector v) {
		if (v.size() < 3) 
			return Matrix();
		if (v.size() >= 3) 
		{
			Matrix A = eye(4, 4);
			for (int i = 0; i < 3; ++i) 
				A(i, A.columns - 1) = v[i];
			return A;
		}
		return Matrix();
	}

	Matrix companionMatrix(Polynomial p) 
	{//Froebenius form
		p = p.makeMonic();
		int N = p.size() - 1;
		std::vector<real> elem;
		elem.resize(N * N, 0);
		int count = 0;
		for (int i = 0; i < N; ++i) 
		{
			for (int j = 0; j < N; ++j) 
			{
				if (j != i - 1 && j < N - 1) 
					elem[i * N + j] = 0;
				if (j == i - 1) 
					elem[i * N + j] = 1;
				if (j == N - 1) 
				{ 
					elem[i * N + j] = -p[count]; 
					count++; 
				}
			}
		}
		return Matrix(N, N, elem);
	}

	Matrix Vandermonde(Polynomial p) 
	{
		// Create a Vandermonde matrix from the vector of coefficients from a_n-1 to a0,
		// which is suitable for finding the complex eigenvalues of the polynomial.
		Matrix CM = companionMatrix(p);
		Vector co(p);
		Matrix V((p.size() - 1), p.size() - 1);
		for (int i = 0; i < p.size() - 1; ++i) 
		{
			for (int j = 0; j < p.size() - 1; ++j) 
			{
				if (!i) 
					V(i, j) = 1;
				else
					V(i, j) = pow(co[i], j + 1);
			}
		}
		return V;
	}

	Matrix outerProduct(Matrix& A, Matrix& B)
	{// Operation that creates for every element of Matrix A, a submatrix of B * (A_i_j)
		if (A.columns != B.rows)
			return Matrix();
		Vector ans(A.size() * B.size());
		int rw = A.rows * B.rows;
		int cl = A.columns * B.columns;

		for (int i = 0; i < rw; ++i) {
			for (int j = 0; j < cl; ++j) {
				ans[i * cl + j] = B.element[((i % B.rows) * B.columns) +
					(j % B.columns)] * A.element[std::floor(i / B.rows) *
					A.columns + std::floor(j / B.columns)];
			}
		}
		return Matrix(rw, cl, ans);
	}

	Matrix wedgeProduct(Matrix& A, Matrix& B)
	{//aka exterior product, ref: https://math.wikia.org/wiki/Outer_product
		Matrix prodA = outerProduct(A, B);
		Matrix prodB = outerProduct(B, A);
		return prodA - prodB;
	}

	Matrix Matrix::multiply(Matrix& A, Matrix& B)
	{
		// Handle case where row/column mismatch prevents multiplication.
		if (A.columns != B.rows)
			return Matrix();

		int size = A.rows * B.columns;
		Vector answerN(size);
		if (B.columns == 1)
		{
			for (int i = 0; i < A.rows; ++i) 
				for (int j = 0; j < A.columns; ++j)
					answerN[i] += A(i, j) * B(j, 0);							

			return Matrix(A.rows, B.columns, answerN);
		}

		for (int i = 0; i < A.rows; ++i)		
			for (int j = 0; j < B.columns; ++j)			
				for (int k = 0; k < A.columns; ++k)				
					answerN[i * B.columns + j] += A(i, k) * B(k, j);								

		return Matrix(A.rows, B.columns, answerN);
	}

	Matrix Matrix::multiply(Matrix& A, real B) {
		Vector answerN(size());
		for (int i = 0; i < A.rows; ++i)
		{
			for (int j = 0; j < A.columns; ++j)
			{
				answerN[(i * A.columns) + j] = A(i, j) * B;
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::multiply(real b, Matrix& A)
	{
		Vector answerN(size());
		for (int i = 0; i < A.rows; ++i)
		{
			for (int j = 0; j < A.columns; ++j)
			{
				answerN[(i * A.columns) + j] = A(i, j) * b;
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::add(Matrix& A, Matrix& B)
	{
		if (A.size() != B.size())
			return Matrix();

		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size, 0);
		for (int i = 0; i < A.rows; ++i)
		{
			for (int j = 0; j < B.columns; ++j)
			{
				answerN[(i * A.columns) + j] = A(i, j) + B(i, j);
			}
		}
		return Matrix(B.rows, A.columns, answerN);
	}

	Matrix Matrix::add(Matrix& A, real B)
	{
		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size, 0);
		for (int i = 0; i < A.rows; ++i)
		{
			for (int j = 0; j < A.columns; ++j)
			{
				answerN[(i * A.columns) + j] = A(i, j) + B;
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::add(real B, Matrix& A)
	{
		return add(A, B);
	}

	Matrix Matrix::subtract(Matrix& A, Matrix& B)
	{
		if (A.size() != B.size())
			return Matrix();
		int size = A.rows * A.columns;
		Vector answerN(size);

		for (int i = 0; i < A.rows; ++i)		
			for (int j = 0; j < B.columns; ++j)			
				answerN[(i * A.columns) + j] += A(i, j) - B(i, j);					

		return Matrix(B.rows, A.columns, answerN);
	}

	Matrix Matrix::subtract(Matrix& A, real B)
	{
		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size, 0);
		for (int i = 0; i < A.rows; ++i)
		{
			for (int j = 0; j < A.columns; ++j)
			{
				answerN[(i * A.columns) + j] = A(i, j) - B;
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::subtract(real b, Matrix& a)
	{
		Matrix temp = -1 * a;
		return add(temp, b);
	}

	Matrix Matrix::divide(Matrix& p, Matrix& q)
	{
		Matrix inv = q.inverse();
		return p * inv;
	}

	Matrix Matrix::divide(Matrix& p2, real q)
	{
		return p2 * (1.0/q);
	}

	Matrix Matrix::divide(real q, Matrix& p2)
	{
		return p2 * (1.0 / q);
	}

	//Overloaded operators:
	bool Matrix::operator==(Matrix& rhs) {
		return (element == rhs.element &&
			rows == rhs.rows && columns == rhs.columns);
	}
	bool Matrix::operator!=(Matrix& rhs) { return !(*this == rhs); }
	Matrix Matrix::operator*=(Matrix& rhs) { *this = *this * rhs; return*this; }
	Matrix Matrix::operator*=(real x) { *this = *this * x; return *this; }
	Matrix Matrix::operator+=(Matrix& rhs) { *this = add(*this, rhs); return*this; }
	Matrix Matrix::operator+=(real x) { *this = add(*this, x); return *this; }
	Matrix Matrix::operator+(Matrix& rhs) { return add(*this, rhs); }
	Matrix Matrix::operator+(real x) { return add(*this, x); }
	Matrix Matrix::operator-=(Matrix& rhs) { *this = subtract(*this, rhs); return*this; }
	Matrix Matrix::operator-=(real x) { *this = subtract(*this, x); return *this; }
	Matrix operator+(Matrix& lhs, Matrix& rhs) { return lhs.add(lhs, rhs); }
	Matrix operator+(Matrix& lhs, real x) { return lhs.add(lhs, x); }
	Matrix operator+(real x, Matrix& rhs) { return rhs.add(x, rhs); }
	Matrix operator-(Matrix& lhs, Matrix& rhs) { return lhs.subtract(lhs, rhs); }
	Matrix operator-(Matrix& lhs, real x) { return lhs.subtract(lhs, x); }
	Matrix operator-(real x, Matrix& rhs) { return rhs.subtract(x, rhs); }
	Matrix operator*(Matrix& lhs, Matrix& rhs) { return lhs.multiply(lhs, rhs); }
	Matrix operator*(Matrix& lhs, real x) { return lhs.multiply(lhs, x); }
	Matrix operator*(real x, Matrix& rhs) { return rhs.multiply(x, rhs); }
	Matrix operator/(Matrix& lhs, Matrix& rhs) { return lhs.divide(lhs, rhs); }
	Matrix operator/(Matrix& lhs, real x) { return lhs.divide(lhs, x); }
	Matrix operator/(real x, Matrix& rhs) { return rhs.divide(x, rhs); }
}