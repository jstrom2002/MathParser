﻿#include "Matrix.h"
#include "Polynomial.h"
#include "Function.h"
#include "Vector.h"
#include "ComplexNumber.h"
#include "MathLib.h"
#include "StringUtils.h"
#include <ctime>
#include <iomanip>

namespace std
{
	MathParser::Matrix pow(MathParser::Matrix m, MathParser::real x) 
	{ 
		return m.exponent(m, x); 
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

	Matrix::Matrix(int row, int col, Vector element2) 
	{
		int n = row;
		int m = col;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row * col, 0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				real n2 = element2[i * col + j];
				element[i * col + j] = n2;
			}
		}
		rows = row; 
		columns = col;
	}

	template <class T>
	Matrix::Matrix(int row, int col, std::vector<T> element2) 
	{
		int n = row;
		int m = col;
		if (col > row) { n = m; }
		if (col < row) { m = n; }
		element.clear();
		element.resize(row * col, 0);
		for (int i = 0; i < row; ++i) {
			for (int j = 0; j < col; ++j) {
				real n2 = element2[i * col + j];
				element[i * col + j] = n2;
			}
		}
		rows = row; 
		columns = col;
	}
	template Matrix::Matrix(int row, int col, std::vector<real> element2);
	template Matrix::Matrix(int row, int col, std::vector<unsigned char> element2);

	Matrix::Matrix(std::vector<Vector> elm) 
	{
		int n = elm.size();
		int m = elm[0].size();
		element.clear();
		element.resize(elm.size() * elm[0].size(), 0);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				real n2 = elm[i][j];
				element[i * m + j] = n2;
			}
		}
		rows = n; columns = m;
	}

	real Matrix::get(int i, int j) 
	{ 
		return element[(i * columns) + j]; 
	}

	void Matrix::set(int i, int j, real x) 
	{ 
		element[(i * columns) + j] = x; 
	}

	int Matrix::size() 
	{ 
		return rows * columns; 
	}

	real Matrix::maxValue() 
	{
		real x = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (std::abs(get(i, j)) > x) { x = get(i, j); }
			}
		}
		return x;
	}

	void Matrix::identity() 
	{
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i == j) 
					element[i * columns + j] = 1;
				else 
					element[i * columns + j] = 0;
			}
		}
	}

	Matrix Matrix::submatrix(int i, int j, int i2, int j2) 
	{
		int size1 = (((i2 + 1) * (j2 + 1)) - ((i + 1) * (j + 1)));
		Vector vals;
		vals.resize(size1,0);
		int counter = 0;
		for (int a = i; a <= i2; ++a) {
			for (int b = j; b <= j2; ++b) {
				vals[counter] = get(a, b);
				++counter;
			}
		}
		return Matrix(i2 - i + 1, j2 - j + 1, vals);
	}

	Vector Matrix::diagonal() 
	{
		Vector vals;
		for (int i = 0; i < rows; ++i) {
			vals.push_back(get(i, i));
		}
		return vals;
	}

	Matrix Matrix::multiply(Matrix A, Matrix B) 
	{
		int size = A.rows * B.columns;
		Vector answerN(size);
		if (B.columns == 1) 
		{
			for (int i = 0; i < A.rows; ++i) {
				for (int j = 0; j < A.columns; ++j) {
					answerN[i] += A.get(i, j) * B.get(j, 0);
				}
			}
			return Matrix(A.rows, B.columns, answerN);
		}

		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < B.columns; ++j) 
			{
				for (int k = 0; k < A.columns; ++k) 
				{
					answerN[(i * B.columns) + j] += A.get(i, k) * B.get(k, j);
				}
			}
		}

		return Matrix(A.rows, B.columns, answerN);
	}

	Matrix Matrix::multiply(Matrix A, real B) {
		Vector answerN(size());
		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < A.columns; ++j) 
			{
				answerN[(i * A.columns) + j] = A.get(i, j) * B; 				
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::add(Matrix A, Matrix B) 
	{
		if (A.size() != B.size())
			return Matrix();

		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size,0);
		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < B.columns; ++j) 
			{
				answerN[(i * A.columns) + j] = A.get(i, j) + B.get(i, j);
			}
		}
		return Matrix(B.rows, A.columns, answerN);
	}

	Matrix Matrix::add(Matrix A, real B) 
	{
		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size,0);
		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < A.columns; ++j) 
			{				
					answerN[(i * A.columns) + j] = A.get(i, j) + B;				
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::subtract(Matrix A, Matrix B) 
	{
		if (A.size() != B.size())
			return Matrix();
		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size,0);
		for (int i = 0; i < size; ++i) { answerN[i] = 0; }

		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < B.columns; ++j) 
			{
				answerN[(i * A.columns) + j] += A.get(i, j) - B.get(i, j);
			}
		}
		return Matrix(B.rows, A.columns, answerN);
	}

	Matrix Matrix::subtract(Matrix A, real B) 
	{
		int size = A.rows * A.columns;
		Vector answerN;
		answerN.resize(size,0);
		for (int i = 0; i < A.rows; ++i) 
		{
			for (int j = 0; j < A.columns; ++j) 
			{
				answerN[(i * A.columns) + j] = A.get(i, j) - B;				
			}
		}
		return Matrix(A.rows, A.columns, answerN);
	}

	Matrix Matrix::inverseExact() 
	{
		if (rows != columns)
			return Matrix();
		Matrix m(rows, columns, element);
		m = m.concatenate(m, identityMatrix(rows, columns));
		m = m.GaussianElimination();
		while (m.columns > columns) { m.removeColumn(0); }
		return m;
	}


	Matrix Matrix::inverse() {
		/*INPUT: n×m matrix A.
		OUTPUT: n×m matrix in reduced row echelon form.
		1.Set j←1
		2.For each row i from 1 to n do
		a. While column j h as all zero elements, set j <- j+1. If j>m return A.
		b. If element a_ij is zero, then interchange row i with a row x > i that has a_xj!=0.
		c. Divide each element of row i by a_ij, thus making the pivot aij equal to one.
		d. For each row k from 1 to n, with k != i, subtract row i multiplied by a_kj from row k.
		3. Return transformed matrix A.		*/

		if (det() == 0) 
			return Matrix();
		Matrix A(rows, columns, element);
		Matrix Aug(rows, columns);
		Aug.identity();

		while (A.zeroInDiag() == true) {
			for (int i = 0; i < rows; ++i) {
				if (A.get(i, i) == 0) {
					for (int j = 0; j < columns; ++j) {
						if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create a non-zero diagonal
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
			real pivot = A.get(j, j);
			for (int i = 0; i < rows; ++i) {
				if (i != j && A.get(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A.get(i, j)) == abs(pivot)) {
						if (A.get(i, j) == pivot) 
						{
							temp = temp * -1.0;
							temp2 = temp2 * -1.0;
							goto skip;
						}
						if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) 
							goto skip;
					}
					if (abs(A.get(i, j)) != abs(pivot)) {
						real a = -A.get(i, j) / pivot;
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
			real pivot = A.get(j, j);
			for (int i = rows - 1; i >= 0; --i) {
				if (i != j && A.get(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A.get(i, j)) == abs(pivot)) {
						if (A.get(i, j) == pivot) {
							temp = temp * -1.0;
							temp2 = temp2 * -1.0;
							goto skip2;
						}
						if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip2; }
					}
					if (abs(A.get(i, j)) != abs(pivot)) {
						real a = -A.get(i, j) / pivot;
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
			if (A.get(i, i) != 0 && A.get(i, i) != 1) {
				real a = ((real)1 / (A.get(i, i)));
				A.multiplyRow(i, a);
				Aug.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, Aug.element);
	}

	Matrix Matrix::exponent(Matrix A, int b) 
	{
		if (b != floor(b)) 
			return Matrix();

		if (b < 0) {
			A = A.inverse();
			b = abs(b);
		}
		if (b == 0) 
		{
			Matrix A2(A.rows, A.columns);
			A2.identity();
			return A2;
		}

		Matrix A2 = A;
		for (int i = 1; i < b; ++i) 
		{ 
			A2 = A2 * A; 
		}

		return A2;
	}

	//Overloaded operators:
	Matrix Matrix::operator*=(Matrix rhs) { *this = multiply(*this, rhs); return*this; }
	Matrix Matrix::operator*=(real x) { *this = multiply(*this, x); return *this; }
	Matrix Matrix::operator*(Matrix rhs) { return multiply(*this, rhs); }
	Matrix Matrix::operator*(real x) { return multiply(*this, x); }

	Matrix Matrix::operator+=(Matrix rhs) { *this = add(*this, rhs); return*this; }
	Matrix Matrix::operator+=(real x) { *this = add(*this, x); return *this; }
	Matrix Matrix::operator+(Matrix rhs) { return add(*this, rhs); }
	Matrix Matrix::operator+(real x) { return add(*this, x); }

	Matrix Matrix::operator-=(Matrix rhs) { *this = subtract(*this, rhs); return*this; }
	Matrix Matrix::operator-=(real x) { *this = subtract(*this, x); return *this; }
	Matrix Matrix::operator-(Matrix rhs) { return subtract(*this, rhs); }
	Matrix Matrix::operator-(real x) { return subtract(*this, x); }


	bool Matrix::operator==(Matrix B) 
	{
		if (rows != B.rows || columns != B.columns) { return false; }
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (get(i, j) != B.get(i, j)) { return false; }
			}
		}
		return true;
	}
	

	Matrix Matrix::outerProduct(Matrix A, Matrix B) {//tensore product that creates for every element of Matrix A, a submatrix of B * (A_i_j)
		if (A.columns != B.rows) { return Matrix(); }
		Vector ans(A.size() * B.size());
		int rw = A.rows * B.rows;
		int cl = A.columns * B.columns;

		for (int i = 0; i < rw; ++i) {
			for (int j = 0; j < cl; ++j) {
				ans[i * cl + j] = B.element[((i % B.rows) * B.columns) + (j % B.columns)] * A.element[floor(i / B.rows) * A.columns + floor(j / B.columns)];
			}
		}
		return Matrix(rw, cl, ans);
	}

	Matrix Matrix::HadamardProduct(Matrix A, Matrix B) {//multiplys two matrices by index so that C(i,j) = A(i,j)*B(i,j)
		if (A.rows != B.rows || A.columns != B.columns) { return Matrix(); }
		Vector elm(A.size());
		for (int i = 0; i < A.rows; ++i) {
			for (int j = 0; j < A.columns; ++j) {
				elm[i * A.columns + j] = A.get(i, j) * B.get(i, j);
			}
		}
		return Matrix(A.rows, A.columns, elm);
	}

	Vector Matrix::row(int rw) 
	{
		if (rw > rows || rw < 0) 		
			return Vector(); 
		
		Vector x;
		x.resize(columns,0);
		for (int k = 0; k < columns; ++k)
		{
			x[k] = get(rw, k);
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
			x[k] = get(k, cl);
		}
		return x;
	}

	void Matrix::setRowNonZeroValues(int rw, real x) 
	{
		for (int k = 0; k < columns; ++k) {
			if (abs(get(rw, k)) > 0) {
				set(rw, k, x);
			}
		}
	}

	void Matrix::setRow(int rw, Vector n) 
	{
		for (int k = 0; k < columns; ++k) 
			set(rw, k, n[k]);
	}

	void Matrix::setColumnNonZeroValues(int col, real x) 
	{
		for (int k = 0; k < rows; ++k) {
			if (abs(get(k, col)) > 0) {
				set(k, col, x);
			}
		}
	}

	void Matrix::setColumn(int col, Vector n) {
		for (int k = 0; k < rows; ++k) { set(k, col, n[k]); }
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
				if (x[j] != get(i, j)) { isRow = false; j = columns; }
				if (j = columns - 1) { isRow = true; index = i; }
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
				if (x[j] != get(i, j)) { isColumn = false; i = rows; }
				if (i = rows - 1) { isColumn = true; index = j; }
			}
		}
		if (index > -1) { return column(index); }
		else { Vector n; return n; }
	}

	void Matrix::reverseRow(int n) 
	{
		Vector vec = row(n);
		std::reverse(vec.begin(), vec.end());
		for (int j = 0; j < columns; ++j) {
			set(n, j, vec[j]);
		}
	}

	void Matrix::reverseColumn(int n) 
	{
		Vector vec = column(n);
		std::reverse(vec.begin(), vec.end());
		for (int j = 0; j < rows; ++j) {
			set(j, n, vec[j]);
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

	Matrix Matrix::expandToUpperLeft(int r, int c) {//takes original matrix and creates a new matrix with the original in the lower right corner.
		Matrix newEl(r, c);
		if (r > rows && c > columns) {
			newEl.identity();
			for (int i = r - rows; i < r; ++i) {
				for (int j = c - columns; j < c; ++j) {
					newEl.set(i, j, get((i - (r - rows)), (j - (c - columns))));
				}
			}
		}
		return Matrix(r, c, newEl.element);
	}

	Matrix Matrix::expandToLowerRight(int r, int c) {//takes original matrix and creates a new matrix with the original in the upper left corner.
		Matrix newEl(r, c);
		if (r > rows && c > columns) {
			newEl.identity();
			for (int i = 0; i < rows; ++i) {
				for (int j = 0; j < columns; ++j) {
					newEl.set(i, j, get(i, j));
				}
			}
		}
		return Matrix(r, c, newEl.element);
	}

	Matrix Matrix::extendRows(int r) {
		Matrix newEl(rows + r, columns);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				newEl.set(i, j, get(i, j));
			}
		}
		return Matrix(rows + r, columns, newEl.element);
	}

	Matrix Matrix::extendColumns(int c) {
		Matrix newEl(rows, columns + c);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				newEl.set(i, j, get(i, j));
			}
		}
		return Matrix(rows, columns + c, newEl.element);
	}

	Matrix Matrix::concatenate(Matrix A, Matrix B) {//creates a larger matrix from two matrices, where each matrix is next to the other one
		int r = A.rows;
		if (B.rows > A.rows) { r = B.rows; }
		int c = A.columns + B.columns;
		Matrix C(r, c);
		for (int i = 0; i < r; ++i) {
			for (int j = 0; j < c; ++j) {
				if (j < A.columns && i < A.rows) { C.set(i, j, A.get(i, j)); }
				if (j >= A.columns && i < B.rows) { C.set(i, j, B.get(i, j - A.columns)); }
			}
		}

		return C;
	}

	Matrix Matrix::autoCorrelation(Matrix A) { return A.convolve(A, A); }//auto-correlation is the convolution of a signal with its own conjugate, which is the function itself in a Real valued domain

	Matrix Matrix::crossCorrelation(Matrix A, Matrix B) {
		int kCenterX = A.columns / 2;
		int kCenterY = A.rows / 2;

		int mm, nn, ii, jj;

		Vector vals(B.element.size());

		for (int i = 0; i < B.rows; ++i) {
			for (int j = 0; j < B.columns; ++j) {
				for (int m = 0; m < A.rows; ++m) {
					for (int n = 0; n < A.columns; ++n) {

						// index of input signal, used for checking boundary
						ii = i + (kCenterY - m);
						jj = j + (kCenterX - n);

						// ignore input samples which are out of bound
						if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) { vals[i * B.columns + j] += B.get(ii, jj) * A.get(m, n); }
					}
				}
			}
		}
		return Matrix(B.rows, B.columns, vals);
	}


	Matrix Matrix::convolve(Matrix A, Matrix B) {
		//covolve Matrix kernel A w/ Matrix B
		int kCenterX = A.columns / 2;
		int kCenterY = A.rows / 2;

		int mm, nn, ii, jj;

		Vector vals(B.element.size());

		for (int i = 0; i < B.rows; ++i) {
			for (int j = 0; j < B.columns; ++j) {
				for (int m = 0; m < A.rows; ++m) {
					mm = A.rows - 1 - m;      // row index of flipped kernel
					for (int n = 0; n < A.columns; ++n) {
						nn = A.columns - 1 - n;  // column index of flipped kernel

						// index of input signal, used for checking boundary
						ii = i + (kCenterY - mm);
						jj = j + (kCenterX - nn);

						// ignore input samples which are out of bound
						if (ii >= 0 && ii < B.rows && jj >= 0 && jj < B.columns) { vals[i * B.columns + j] += B.get(ii, jj) * A.get(mm, nn); }
					}
				}
			}
		}
		return Matrix(B.rows, B.columns, vals);
	}

	bool Matrix::isIntegerMatrix() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (floor(abs(get(i, j))) != abs(get(i, j)))
					return false;
			}
		}
		return true;
	}

	bool Matrix::isSymmetric() {
		if (*this == transpose()) 
			return true;
		return false;
	}

	bool Matrix::isBoolean() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (get(i, j) != 0 && get(i, j) != 1) 
					return false;
			}
		}
		return true;
	}

	bool Matrix::isDiagonalized() {
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i != j) {
					if (abs(get(i, j)) > 0) 
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
		if (det() == 0) 
			return false;		
		return true;
	}

	bool Matrix::isLinearlyIndependent() {
		if (rows < columns) { return true; }
		if (det() == 0) { return false; }
		if (GramMatrix().det() == 0) { return false; }
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
				element[i*rows + j] = get(j, i);
			}
		}
		return Matrix();
	}

	Matrix Matrix::populationCovarianceMatrix() 
	{
		Matrix M(columns, columns);
		for (int i = 0; i < columns; ++i) {
			for (int j = 0; j < columns; ++j) {
				Vector X1 = column(i) - columnMean(i);
				Vector X2 = column(j) - columnMean(j);
				M.set(i, j, dot(X1, X2) / (rows));
			}
		}
		return M;
	}

	Matrix Matrix::sampleCovarianceMatrix() 
	{
		Matrix M(columns, columns);
		for (int i = 0; i < columns; ++i) {
			for (int j = 0; j < columns; ++j) {
				Vector X1 = column(i) - columnMean(i);
				Vector X2 = column(j) - columnMean(j);
				M.set(i, j, dot(X1, X2) / ((rows)-1));
			}
		}
		return M;
	}


	Matrix Matrix::transpose() {
		int n = rows;
		int m = columns;
		Vector p(n * m);
		for (int i = 0; i < m; ++i) {
			for (int j = 0; j < n; ++j) {
				p[i * rows + j] = get(j, i);
			}
		}
		return Matrix(m, n, p);
	}

	Matrix Matrix::GramMatrix() {//turn row vector or nxn matrix into its Gram matrix
		Matrix T = transpose();
		if (columns == 1) { return multiply(*this, T); }
		else { return multiply(T, *this); }
	}

	Matrix Matrix::GramMatrix(Matrix M) {//turn 2 vectors into a Gram matrix
		Matrix Trans = M.transpose();
		if (M.columns == 1) { return multiply(M, Trans); }
		else { return multiply(Trans, M); }
	}

	Matrix Matrix::Vandermonde() {//create Vandermonde matrix from vector
		int p = columns;
		if (columns == 1) { p = rows; }
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
				answer += get(i, j);
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
				answer *= get(i, j);
			}
		}
		return pow(answer, real(1.0 / size()));
	}

	real Matrix::trace() 
	{
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += get(i, i);
		}
		return answer;
	}

	std::string Matrix::to_string(int precision) {
		std::ostringstream s;
		int defaultLength = 4;
		int largestElement = this->largestElement();
		bool tf = isIntegerMatrix();
		if (tf == true) { --defaultLength; }
		if (tf == false) { ++defaultLength; }
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				s <<
					std::fixed <<
					std::left <<
					std::setw(defaultLength + (precision)+(largestElement)) <<
					std::setprecision(precision) <<
					element[i * columns + j];
			}
			s << L"\r\n";
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
				answer += pow(get(i, j), 2);
			}
		}
		return sqrt(answer);
	}

	int Matrix::largestElement()
	{
		int lelem = 0;
		for (int i = 0; i < rows; ++i) 
		{
			for (int j = 0; j < columns; ++j) 
			{
				int len = getLength(element[i*columns+j]);
				if (len > lelem)
					lelem = len; 				
			}
		}
		return lelem;
	}

	real Matrix::FrobeniusNorm() 
	{ //the Froebenius norm is equal to the sqrt of the trace of 
	  //the matrix multiplied by its conjugate transpose.
		Matrix AT = *this;
		AT.transpose();
		return sqrt(multiply(*this, AT).trace());
	}

	real Matrix::pNorm(real p) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(get(i, j), p);
			}
		}
		return sqrt(answer);
	}

	real Matrix::sumRow(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += get(r, i);
		}
		return answer;
	}

	real Matrix::sumColumn(int r) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += get(i, r);
		}
		return answer;
	}

	real Matrix::sumAll() {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; j++) {
				answer += get(i, j);
			}
		}
		return answer;
	}

	real Matrix::columnMean(int c) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += get(i, c);
		}
		return (answer / rows);
	}

	real Matrix::rowMean(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += get(r, i);
		}
		return (answer / columns);
	}

	real Matrix::rowNorm(int r) {
		real answer = 0;
		for (int i = 0; i < columns; ++i) {
			answer += pow(get(r, i), 2);
		}
		return sqrt(answer);
	}

	real Matrix::columnNorm(int c) {
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(i, c), 2);
		}
		return sqrt(answer);
	}

	real Matrix::columnSumSquares(int c) {
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(i, c) - colMean, 2);
		}
		return answer;
	}

	real Matrix::rowSumSquares(int r) {
		real rMean = rowMean(r);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(r, i) - rMean, 2);
		}
		return answer;
	}

	real Matrix::columnVariance(int c) {
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(i, c) - colMean, 2);
		}
		return answer / (rows - 1);
	}

	real Matrix::columnCovariance(int c1, int c2) {
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

	real Matrix::columnDeviation(int c) {
		real colMean = columnMean(c);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(i, c) - colMean, 2);
		}
		return sqrt(answer / (rows - 1));
	}

	real Matrix::rowCovariance(int r1, int r2) {
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
			answer += pow(get(r, i) - rMean, 2);
		}
		return answer / (columns - 1);
	}

	real Matrix::rowDeviation(int r) {
		real rMean = rowMean(r);
		real answer = 0;
		for (int i = 0; i < rows; ++i) {
			answer += pow(get(r, i) - rMean, 2);
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
				answer += pow(get(i, j) - Mean, 2);
			}
		}
		return sqrt(answer / ((rows * columns) - 1));
	}

	real Matrix::populationStandardDeviation() {
		real answer = 0;
		real Mean = mean();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				answer += pow(get(i, j) - Mean, 2);
			}
		}
		return sqrt(answer / (rows * columns));
	}

	real Matrix::ChiSquareTestStatistic() {
		real answer = 0;
		real df = (rows - 1) * (columns - 1);

		real sm = sumAll();
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real Ei = (sumRow(i) * sumColumn(j) / sm);
				answer += pow(get(i, j) - Ei, 2) / Ei;
			}
		}
		return answer;
	}

	real Matrix::ChiSquareDegreesFreedom() { return  (rows - 1) * (columns - 1); }

	void Matrix::removeRow(int n) {
		Vector newElements;
		newElements.resize((rows - 1) * columns, 0);
		int counter = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i != n) { newElements[counter] = get(i, j); ++counter; }
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
				if (j != n) { newElements[counter] = get(i, j); ++counter; }
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
		Matrix T(A.rows + B.rows, A.columns + B.columns);
		T.identity();
		for (int i = 0; i < T.rows; ++i) {
			for (int j = 0; j < T.columns; ++j) {
				if (i < A.rows && j < A.columns) { T.set(i, j, A.get(i, j)); }
				if (i >= A.rows && j >= A.columns) { T.set(i, j, B.get(i - A.rows, j - A.columns)); }
			}
		}
		return T;
	}

	Matrix Matrix::tensorProduct(Matrix A, Matrix B) {
		Matrix T(A.rows * B.rows, A.columns * B.columns);
		T.identity();
		int Rbound = (A.rows - 1);
		int Cbound = (A.columns - 1);
		int Rcounter = 0;
		int Ccounter = 0;
		for (int i = 0; i < T.rows; ++i) {
			for (int j = 0; j < T.columns; ++j) {
				T.set(i, j, A.get(Rcounter, Ccounter) * B.get(i % B.rows, j % B.columns));
				if (j % A.columns == Cbound && j != 0) { ++Ccounter; }
			}
			Ccounter = 0;
			if (i % A.rows == Rbound && i != 0) { ++Rcounter; }
		}
		return T;
	}

	bool Matrix::canSwapOutZeroDiagonals() {
		for (int i = 0; i < columns; ++i) {
			if (sumColumn(i) == 0) { return false; }
		}
		return true;
	}

	Matrix Matrix::swapOutZeroDiagonals() {
		Matrix M(rows, columns, element);
		while (M.zeroInDiag() == true) {
			for (int i = 0; i < rows; ++i) {
				if (M.get(i, i) == 0) {
					for (int j = 0; j < columns; ++j) {
						if (M.get(i, j) != 0 && M.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create M non-zero diagonal
							M.swapRows(i, j); j = columns;
						}
					}
				}
			}
		}
		return M;
	}

	real Matrix::columnAbsMax(int n) {
		real m = get(0, n);
		if (rows <= 1) { return m; }
		for (int i = 0; i < rows; ++i) {
			real temp = (get(i, n));
			if (abs(temp) > m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnMax(int n) {
		real m = get(0, n);
		if (rows <= 1) { return m; }
		for (int i = 0; i < rows; ++i) {
			real temp = get(i, n);
			if (temp > m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnAbsMin(int n) {
		real m = get(0, n);
		if (rows <= 1) { return m; }
		for (int i = 1; i < rows; ++i) {
			real temp = (get(i, n));
			if (abs(temp) < m) { m = temp; }
		}
		return m;
	}

	real Matrix::columnMin(int n) {
		real m = get(0, n);
		if (rows <= 1) { return m; }
		for (int i = 1; i < rows; ++i) {
			real temp = get(i, n);
			if (temp < m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowAbsMax(int n) {
		real m = get(n, 0);
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = (get(n, i));
			if (abs(temp) > m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowMax(int n) {
		real m = get(n, 0);
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = get(n, i);
			if (temp > m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowAbsMin(int n) {
		real m = get(n, 0);
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = (get(n, i));
			if (abs(temp) < m) { m = temp; }
		}
		return m;
	}

	real Matrix::rowMin(int n) {
		real m = get(n, 0);
		if (columns <= 1) { return m; }
		for (int i = 0; i < columns; ++i) {
			real temp = get(n, i);
			if (temp < m) { m = temp; }
		}
		return m;
	}

	int Matrix::rowNonzeroValues(int rw) {
		int n = 0;
		for (int i = 0; i < columns; ++i) {
			if (abs(get(rw, i)) > 0) {
				++n;
			}
		}
		return n;
	}

	int Matrix::columnNonzeroValues(int cl) {
		int n = 0;
		for (int i = 0; i < rows; ++i) {
			if (abs(get(i, cl)) > 0) {
				++n;
			}
		}
		return n;
	}

	int Matrix::getPivot(int rw) {
		int val = 0;
		while (val < columns) {
			if (abs(get(rw, val)) > 0) {
				return val;
			}
			++val;
		}
		return -1;
	}

	int Matrix::getReversePivot(int rw) {
		int val = columns - 1;
		while (val >= 0) {
			if (abs(get(rw, val)) > 0) {
				return val;
			}
			--val;
		}
		return -1;
	}

	Matrix Matrix::GaussianElimination() {
		int i, j, k;
		real piv;
		Matrix m(rows, columns, element);
		//Gauss-Jordan elimination
		for (i = 0; i < rows; i++) {
			if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
			int ind = m.getPivot(i);
			if (ind >= 0) {
				piv = m.get(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) {
					for (int l = i; l < rows; ++l) {
						m.set(i, l, m.get(i, l) / piv);
					}
					m.set(i, i, 1);
				}
				//proceed down the column to make each non-pivot value = 0
				for (j = 0; j < rows; j++) {
					if (m.get(j, i) != 0) {
						if (j != i) {
							piv = m.get(j, i);
							for (k = i + 1; k < columns; k++) {
								m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
							}
							m.set(j, i, 0);
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
			if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
			int ind = m.getReversePivot(i);
			if (ind >= 0) {
				piv = m.get(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) {
					for (int l = 0; l < i; ++l) {
						m.set(i, l, m.get(i, l) / piv);
					}
					m.set(i, i, 1);
				}
				//proceed up the column to make each non-pivot value = 0
				for (j = i - 1; j >= 0; --j) {
					if (m.get(j, i) != 0) {
						if (j != i) {
							piv = m.get(j, i);//m.get(i,i);
							for (k = 0; k < i; ++k) {
								m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
							}
							m.set(j, i, 0);
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
			if (zeroInDiag() && canSwapOutZeroDiagonals()) { m = m.swapOutZeroDiagonals(); }//remove zeroes in diagonal if possible
			int ind = m.getPivot(i);
			if (ind >= 0) {
				piv = m.get(i, ind);
				//normalize pivot row so pivot = 1
				if (piv != 1 && piv != 0 && abs(piv) > 0) {
					for (int l = i + 1; l < rows; ++l) {
						m.set(i, l, m.get(i, l) / piv);
					}
					m.set(i, i, 1);
				}
				//proceed down the column to make each non-pivot value = 0
				for (j = i + 1; j < rows; j++) {
					if (m.get(j, i) != 0) {
						if (j != i) {
							piv = m.get(j, i);//m.get(i,i);
							for (k = i + 1; k < columns; k++) {
								m.set(j, k, m.get(j, k) - (m.get(i, k) * piv));
							}
							m.set(j, i, 0);
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

	int Matrix::nullity() { return columns - rank(); }

	void Matrix::removeZeroColumns() {
		for (int i = 0; i < columns; ++i) {
			if (sumColumn(i) == 0) { removeColumn(i); --i; }
		}
	}

	void Matrix::removeZeroRows() {
		for (int i = 0; i < rows; ++i) {
			if (sumRow(i) == 0) {
				removeRow(i); --i;
			}
		}
	}

	std::vector<int> Matrix::pivotColumns() {
		std::vector<int> piv;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				real val = get(i, j);
				if (val != 0) {
					piv.push_back(j);
					j = columns + 1;
				}
			}
		}
		return piv;
	}

	Matrix Matrix::nullSpace() {
		Matrix RRE = GaussianElimination();
		int kernelDimension = RRE.rank();//the rank of the matrix is equivalent to the dimension of the kernel
		//if (columns == kernelDimension) { return Matrix(1, 1); }

		//for a RRE matrix, the non-pivot columns are the free variables
		std::vector<int> piv = RRE.pivotColumns();
		RRE.removeZeroRows();

		std::vector<Vector> ans;
		int counter = 0;
		while (ans.size() < kernelDimension && counter < columns * 10) {
			for (int i = 0; i < piv.size(); ++i) {
				Matrix Mcopy(RRE);//set up matrix

				//make one of the pivot columns = 1, the others = 0. rotate which one is 1 over each iteration to find all distinct solution vectors
				for (int p = 0; p < piv.size(); ++p) {
					if (p != i) { Mcopy.setColumnNonZeroValues(piv[p], 0); }
				}
				Mcopy.setColumnNonZeroValues(piv[i], 1);

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
			real pivot = M.get(j + 1, j);
			for (int i = j + 2; i < rows; ++i) {
				Vector temp = M.row(j + 1);
				if (M.get(i, j) != 0) {
					if (abs(M.get(i, j)) == abs(pivot)) {
						if (M.get(i, j) == pivot) { temp*=-1; goto skipUH; }
						if ((M.get(i, j) == (-1 * pivot)) || ((-1 * M.get(i, j)) == pivot)) { goto skipUH; }
					}
					if (abs(M.get(i, j)) != abs(pivot)) { temp*=(-M.get(i, j) / pivot); }

				skipUH:
					M.addToRow(i, temp);
				}
			}
		}
		return Matrix(rows, columns, M.element);
	}

	Matrix Matrix::lowerHessenbergForm() {
		if (rows != columns) { return Matrix(); }
		Matrix A(rows, columns, element);

		for (int j = columns - 1; j >= 0; --j) {//go the other way to remove the upper triangle of non-diagonal values
			real pivot = A.get(j, j);
			for (int i = j - 1; i >= 0; --i) {
				if (i != j && A.get(i, j) != 0) {
					Vector temp = A.row(j);

					if (abs(A.get(i, j)) == abs(pivot)) {
						if (A.get(i, j) == pivot) {
							temp*= -1;
							goto skipLH;
						}
						if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) 
						{ 
							goto skipLH; 
						}
					}
					if (abs(A.get(i, j)) != abs(pivot)) {
						real a = -A.get(i, j) / pivot;
						temp*= a;
					}

				skipLH:
					A.addToRow(i, temp);
				}
			}
		}

		return Matrix(rows, columns, A.element);
	}

	Matrix Matrix::characteristicMatrix() {//numerical method
		Matrix next(rows, columns, element);
		for (int k = 1; k < rows; ++k) {
			Matrix M(columns, rows);
			M.identity();
			for (int p = 1; p < columns + 1; ++p) {
				if (p == (rows - k)) {
					M.element[((rows - k - 1) * rows) + (p - 1)] = (1 / (next.element[((rows - k) * columns) + (rows - k - 1)]));
				}
				else { M.element[((rows - k - 1) * rows) + (p - 1)] = -next.element[((rows - k) * columns) + (p - 1)] / (next.element[((columns - k) * columns) + (columns - k - 1)]); }
			}
			Matrix Minv(columns, rows);
			Minv.identity();
			for (int p = 1; p < columns + 1; ++p) { Minv.element[(columns - k - 1) * columns + (p - 1)] = next.element[(columns - k) * columns + (p - 1)]; }
			next = next.multiply(Minv, next.multiply(next, M));
		}

		for (int i = 1; i < next.rows; ++i) {//blank out anything that's not in the first row;
			for (int j = 0; j < next.columns; ++j) {
				if (i == j) { next.element[i * columns + j] = 1; }
				else { next.element[i * columns + j] = 0; }
			}
		}
		return Matrix(rows, columns, next.element);
	}


	bool Matrix::zeroInDiag() {//test if there are any 0's in the diagonal of a matrix
		for (int i = 0; i < rows; ++i) {
			if (get(i, i) == 0) { return true; }
		}
		return false;
	}

	real Matrix::det() {//determinant by exact (ie non-numerical) method, achieved by taking the product of diagonals in the reduced row echelon form (via Gaussian elimination)
		if (rows != columns) { return 0; }
		if (rows == 1) {
			return get(0, 0);
		}

		//MATRICES OF 2nd ORDER
		if (rows == 2) {
			return (get(0, 0) * get(1, 1) - get(0, 1) * get(1, 0));
		}

		//NTH ORDER -- use elimination
		Matrix M(rows, columns, element);

		int signFlip = 0;
		while (M.zeroInDiag() == true) {
			for (int i = 0; i < rows; ++i) {
				if (M.get(i, i) == 0) {
					for (int j = 0; j < columns; ++j) {
						if (M.get(i, j) != 0 && M.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create M non-zero diagonal
							M.swapRows(i, j); j = columns;
							++signFlip;	//each time you swap rows, you flip the sign of the determinant
						}
					}
				}
			}
		}

		//elimination with 3 cases:  
		for (int j = 0; j < columns; ++j) {
			real pivot = M.get(j, j);
			for (int i = j + 1; i < rows; ++i) {
				Vector temp = M.row(j);
				if (M.get(i, j) != 0) {
					if (abs(M.get(i, j)) == abs(pivot)) {
						if (M.get(i, j) == pivot) 
						{ 
							temp*= -1; 
							goto skip; 
						}
						if ((M.get(i, j) == (-1 * pivot)) || ((-1 * M.get(i, j)) == pivot)) { goto skip; }
					}
					if (abs(M.get(i, j)) != abs(pivot)) 
					{ 
						temp*= (-M.get(i, j) / pivot); 
					}

				skip:
					M.addToRow(i, temp);
				}
			}
		}

		real answer = 1;
		for (int i = 0; i < columns; ++i) {
			answer *= M.get(i, i);//multiply diagonals to get determinant
		}
		answer *= pow(-1, signFlip % 2);//fix flipped determinant
		return answer;
	}//END DETERMINANT FUNCTION=======================================

	real Matrix::detExact() {
		Matrix m = upperTriangularize();
		real det = 1;
		for (int i = 0; i < m.rows; ++i) {
			det *= m.get(i, i);
		}
		return det;
	}

	Matrix Matrix::reducedRowEchelonForm() {
		Matrix A(rows, columns, element);

		while (A.zeroInDiag() == true) {
			for (int i = 0; i < rows; ++i) {
				if (A.get(i, i) == 0) {
					for (int j = 0; j < columns; ++j) {
						if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) { //carefully choose rows to swap to create a non-zero diagonal
							A.swapRows(i, j);
							j = columns;
						}
					}
				}
			}
		}

		for (int j = 0; j < rows; ++j) {//go down and eliminate the lower triangle of values
			real pivot = A.get(j, j);
			for (int i = 0; i < rows; ++i) {
				if (i != j && A.get(i, j) != 0) {
					Vector temp = A.row(j);

					if (abs(A.get(i, j)) == abs(pivot)) {
						if (A.get(i, j) == pivot) {
							temp*= -1;
							goto skip;
						}
						if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip; }
					}
					if (abs(A.get(i, j)) != abs(pivot)) {
						real a = -A.get(i, j) / pivot;
						temp*= a;
					}

				skip:
					A.addToRow(i, temp);

				}
			}
		}


		for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
			if (A.get(i, i) != 0 && A.get(i, i) != 1) {
				real a = ((real)1 / (A.get(i, i)));
				A.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, A.element);
	}

	Matrix Matrix::addRow(Vector vec) {//adds vector to right hand side by default
		Matrix M(rows + 1, columns);
		int rw = rows;
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (i == rw) {
					M.set(i, j, vec[j]);
				}
				else {
					if (i < rw) {
						M.set(i, j, get(i, j));
					}
					if (i > rw) {
						M.set(i, j, get(i - 1, j));
					}
				}
			}
		}
		return M;
	}

	Matrix Matrix::addRow(int rw, Vector vec) {
		Matrix M(rows + 1, columns);
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (i == rw) {
					M.set(i, j, vec[j]);
				}
				else {
					if (i < rw) {
						M.set(i, j, get(i, j));
					}
					if (i > rw) {
						M.set(i, j, get(i - 1, j));
					}
				}
			}
		}
		return M;
	}

	Matrix Matrix::addColumn(Vector vec) {//adds to bottom of matrix by default
		Matrix M(rows, columns + 1);
		int cl = columns;
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (j == cl) {
					M.set(i, j, vec[i]);
				}
				else {
					if (j < cl) {
						M.set(i, j, get(i, j));
					}
					if (j > cl) {
						M.set(i, j, get(i, j - 1));
					}
				}
			}
		}
		return M;
	}

	Matrix Matrix::addColumn(int cl, Vector vec) {
		Matrix M(rows, columns + 1);
		for (int i = 0; i < M.rows; ++i) {
			for (int j = 0; j < M.columns; ++j) {
				if (j == cl) {
					M.set(i, j, vec[i]);
				}
				else {
					if (j < cl) {
						M.set(i, j, get(i, j));
					}
					if (j > cl) {
						M.set(i, j, get(i, j - 1));
					}
				}
			}
		}
		return M;
	}

	bool Matrix::isInconsistent(Vector b) {//by Rouche-Capelli Theorem,  any system of linear equations (underdetermined or otherwise) is inconsistent if the rank of the augmented matrix is greater than the rank of the coefficient matrix
		if (rank() <= addColumn(columns - 1, b).rank()) { return true; }
		return false;
	}

	bool Matrix::isUnderDetermined() { return (rows < columns); }
	bool Matrix::isOverDetermined() { return (rows >= columns); }

	bool Matrix::isSingular() {

		//FILL IN
		return false;
	}

	bool Matrix::solutionCheck(Vector x, Vector b) {//tests if Ax=b is a solution
		for (int i = 0; i < rows; ++i) {
			real val = 0;
			for (int j = 0; j < columns; ++j) {
				val += get(i, j) * x[j];
			}
			if (abs(val) - abs(b[i]) > 0) { return false; }
		}
		return true;
	}

	Vector Matrix::solve(Vector b) {//solve for Ax=b, where x will be the returned vector of variable values
		if (isUnderDetermined() == true) {

		}
		if (isOverDetermined() == true) {

		}

		//if m==n
		int n = rows;
		Vector x(n);
		for (int i = n - 1; i >= 0; i--) {
			x[i] = get(i, n - 1) / get(i, i);
			for (int k = i - 1; k >= 0; k--) {
				set(k, n - 1, get(k, n - 1) - get(k, i) * x[i]);
			}
		}
		return x;

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
		Matrix MMT = M.multiply(M, MT);
		Mat.push_back(MMT.eigenvectors());

		//search for V
		Matrix MTM = M.multiply(MT, M);
		Mat.push_back(MTM.eigenvectors());

		//search for S
		eigens1 = MMT.eigenvaluesRealExact();
		eigens2 = MTM.eigenvaluesRealExact();
		eigens1.insert(eigens1.begin(), eigens2.begin(), eigens2.end());
		eigens2.clear();
		eigens1.removeNullValues();
		Mat.push_back(Vector(eigens1).toDiagonalMatrix());
		return(Mat);
	}

	std::vector<Matrix> Matrix::LU() {
		std::vector<Matrix> answer;
		Matrix A(rows, columns, element);
		Matrix L(rows, columns);
		Matrix U(rows, columns);

		int i, j, k;
		real sum = 0;

		for (i = 0; i < rows; i++) {
			U.set(i, i, 1);
		}

		for (j = 0; j < columns; j++) {
			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(i, k) * U.get(k, j);
				}
				L.set(i, j, A.get(i, j) - sum);
			}

			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(j, k) * U.get(k, i);
				}
				if (L.get(j, j) == 0) {
					//printf("det(L) may be close to 0!\n Can't divide by 0...\n");
					answer[0] = Matrix();
					answer[1] = Matrix();
					return answer;
				}
				U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
			}
		}
		answer.push_back(L);
		answer.push_back(U);
		return answer;
	}

	Matrix Matrix::L() {
		Matrix A(rows, columns, element);
		Matrix L(rows, columns);
		Matrix U(rows, columns);

		int i, j, k;
		real sum = 0;

		for (i = 0; i < rows; i++) {
			U.set(i, i, 1);
		}

		for (j = 0; j < columns; j++) {
			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(i, k) * U.get(k, j);
				}
				L.set(i, j, A.get(i, j) - sum);
			}

			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(j, k) * U.get(k, i);
				}
				if (L.get(j, j) == 0) {
					printf("det(L) close to 0!\n Can't divide by 0...\n");
					return Matrix();
				}
				U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
			}
		}
		return L;
	}

	Matrix Matrix::U() {
		Matrix A(rows, columns, element);
		Matrix L(rows, columns);
		Matrix U(rows, columns);

		int i, j, k;
		real sum = 0;

		for (i = 0; i < rows; i++) {
			U.set(i, i, 1);
		}

		for (j = 0; j < columns; j++) {
			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(i, k) * U.get(k, j);
				}
				L.set(i, j, A.get(i, j) - sum);
			}

			for (i = j; i < rows; i++) {
				sum = 0;
				for (k = 0; k < j; k++) {
					sum = sum + L.get(j, k) * U.get(k, i);
				}
				if (L.get(j, j) == 0) {
					printf("det(L) close to 0!\n Can't divide by 0...\n");
					return Matrix();
				}
				U.set(j, i, ((A.get(j, i) - sum) / L.get(j, j)));
			}
		}
		return U;
	}

	std::vector<Matrix> Matrix::QR() {//QR decomposition by Gram-Schmidt process
		std::vector<Matrix> A;
		Matrix Q(rows, columns);
		Matrix R(rows, columns);

		//set up Q matrix
		//===============
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
		//===============
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (j >= i) {
					Vector a = column(j);
					Vector b = Q.column(i);
					real temp = dot(a, b);
					R.set(i, j, temp);
				}
			}
		}
		A.push_back(Q);
		A.push_back(R);
		return A;
	}

	Matrix Matrix::Q() {//QR decomposition by Gram-Schmidt process
		Matrix Q(rows, columns);
		Matrix R(rows, columns);

		//set up Q matrix
		//===============
		Q.setColumn(0, column(0).length());
		for (int i = 1; i < rows - 1; ++i) {
			Vector c = column(i);
			for (int j = i; j > 0; --j) {
				Vector e = Q.column(j - 1);
				real dp = dot(c, e);
				Vector tp = (e* dp);
				c = (c - tp);
			}
			c = c.length();
			Q.setColumn(i, c);
		}

		//set up R matrix
		//===============
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (j >= i) {
					Vector a = column(j);
					Vector b = Q.column(i);
					real temp = dot(a, b);
					R.set(i, j, temp);
				}
			}
		}
		return Matrix(rows, columns, Q.element);
	}

	Matrix Matrix::R() {//QR decomposition by Gram-Schmidt process
		Matrix Q(rows, columns);
		Matrix R(rows, columns);

		//set up Q matrix
		//===============
		Q.setColumn(0, column(0).length());
		for (int i = 1; i < rows - 1; ++i) {
			Vector c = column(i);
			for (int j = i; j > 0; --j) {
				Vector e = Q.column(j - 1);
				real dp = dot(c, e);
				Vector tp = (e* dp);
				c -= tp;
			}
			c = c.length();
			Q.setColumn(i, c);
		}

		//set up R matrix
		//===============
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (j >= i) {
					Vector a = column(j);
					Vector b = Q.column(i);
					real temp = dot(a, b);
					R.set(i, j, temp);
				}
			}
		}
		return Matrix(rows, columns, R.element);
	}

	Matrix Matrix::inverseByQR() {
		if (det() == 0) { return Matrix(); }
		std::vector<Matrix> qr = QR();
		Matrix A = multiply(qr[1], qr[0].transpose());
		return Matrix(rows, columns, A.element);
	}

	Matrix Matrix::GramSchmidt() {//orthogonalization by Gram-Schmidt process using Gaussian Elimination

		if (det() == 0) { return Matrix(); }
		Matrix Aug(rows, columns, element);
		Matrix AT(rows, columns, element);
		AT = AT.transpose();
		Matrix A = multiply(*this, AT);

		while (A.zeroInDiag()) 
		{
			for (int i = 0; i < rows; ++i) 
			{
				if (A.get(i, i) == 0) 
				{
					for (int j = 0; j < columns; ++j) 
					{
						//carefully choose rows to swap to create a non-zero diagonal
						if (A.get(i, j) != 0 && A.get(j, i) != 0 && i != j) 
						{
							A.swapRows(i, j);
							Aug.swapRows(i, j);
							j = columns;
						}
					}
				}
			}
		}

		for (int j = 0; j < columns; ++j) {//go down and eliminate the lower triangle of values
			real pivot = A.get(j, j);
			for (int i = j + 1; i < rows; ++i) {
				if (i != j && A.get(i, j) != 0) {
					Vector temp = A.row(j);
					Vector temp2 = Aug.row(j);

					if (abs(A.get(i, j)) == abs(pivot)) {
						if (A.get(i, j) == pivot) {
							temp *= -1;
							temp2 *= -1;
							goto skip;
						}
						if ((A.get(i, j) == (-1 * pivot)) || ((-1 * A.get(i, j)) == pivot)) { goto skip; }
					}
					if (abs(A.get(i, j)) != abs(pivot)) 
					{
						real a = -A.get(i, j) / pivot;
						temp *= a;
						temp2 *= a;
					}

				skip:
					A.addToRow(i, temp);
					Aug.addToRow(i, temp2);

				}
			}
		}

		for (int i = 0; i < rows; ++i) {//make diagonals into all 1's
			if (A.get(i, i) != 0 && A.get(i, i) != 1) {
				real a = ((real)1 / (A.get(i, i)));
				A.multiplyRow(i, a);
				Aug.multiplyRow(i, a);
			}
		}
		return Matrix(rows, columns, Aug.element);
	}

	std::vector<int> Matrix::JacobiIndexing() {//returns indices i,j for largest element in matrix quickly for use in JacobiTransformation method
		std::vector<int> answer;
		answer.push_back(-1);
		answer.push_back(-1);
		real largest = 0;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i != j && abs(get(i, j)) > largest) {
					largest = abs(get(i, j));
					answer[0] = i;
					answer[1] = j;
				}
			}
		}
		return answer;
	}

	Matrix Matrix::GivensRotationMatrix(real a1, real a2, int r1, int r2) {//eliminate a2 with a1, a1 is in row r1, a2 is in row r2
		Matrix el(rows, rows);
		el.identity();

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

	Matrix Matrix::JacobiRotationMatrix(int p, int q, real c, real s) {//used in Jacobi Transformation -- will diagonalize a matrix and produce an eigenvector matrix
		Matrix el(rows, rows);
		el.identity();

		el.element[p * rows + p] = c;
		el.element[q * rows + q] = c;
		el.element[p * rows + q] = s;
		el.element[q * rows + p] = -s;
		return Matrix(rows, rows, el.element);
	}

	std::vector<Matrix> Matrix::JacobiTransformation() {//outputs an array with both diagonalized matrix and eigenvectors
														//if (rows != columns || isSymmetric()==false) { return NULL;	}
		Matrix D(rows, columns, element);  //matrix to rotate/diagonalize
		Matrix Eigenvectors(rows, columns); //matrix of Eigenvectors produced by multiplying each successive Jacobi Rotation matrix
		Eigenvectors.identity();

		int counter = 0;//counts the number of iterations
		long unsigned int limit = 10000;
		while (D.isDiagonalized() == false && counter < limit) {
			std::vector<int> elem = D.JacobiIndexing();
			if (elem[0] > D.rows || elem[1] > D.columns || elem[0] < 0 || elem[1] < 0 || elem[0] == elem[1]) { break; }
			int p = elem[0];
			int q = elem[1];
			real angle = (D.get(q, q) - D.get(p, p)) / (2 * D.get(p, q));
			int sign = 1;
			if (angle < 0) { sign = -1; }
			real t = sign / (abs(angle) + sqrt((angle * angle) + 1));
			real c = 1 / sqrt((t * t) + 1);
			real s = t * c;
			Matrix rot = D.JacobiRotationMatrix(p, q, c, s);
			Matrix rotT = rot.transpose();
			D = D.multiply(rotT, D);
			D *= rot;
			Eigenvectors = Eigenvectors.multiply(Eigenvectors, rot);
			++counter;
		}
		std::vector<Matrix> answer;
		answer.push_back(D);
		answer.push_back(Eigenvectors);
		return answer;
	}

	Matrix Matrix::eigenvectors() {
		if (isSymmetric() == true) {
			std::vector<Matrix> EV = JacobiTransformation();
			return EV[1];
		}
		else {
			std::vector<Matrix> qr = QR();
			Matrix Eigenvalues = qr[0];
			Matrix A2 = multiply(qr[1], qr[0]);

			int iterations = pow(rows * columns, 3);
			if (iterations > 500) { iterations = 500; }
			for (int i = 1; i < iterations; ++i) {//then, QR decomposition. Use Q as multiplicand where Eigenvectors = Q1*Q2*...Qn
				qr = A2.QR();
				if (qr[0].rows == 0 && qr[1].rows == 0) { break; }
				Eigenvalues = Eigenvalues.multiply(Eigenvalues, qr[0]);
				A2 = A2.multiply(qr[1], qr[0]);
				qr.clear();
			}
			return Eigenvalues;
		}
		return Matrix();
	}

	Matrix Matrix::Cholesky() {	//returns matrix L decomposition of matrix where A = LL^T
		if (rows <= 1 || columns <= 1) { return Matrix(); }
		Matrix L(rows, columns);
		L.set(0, 0, sqrt(get(0, 0)));

		for (int j = 0; j < columns; ++j) {
			for (int i = j; i < rows; ++i) {
				if (i == j) {
					real temp = get(j, j);
					real temp2 = 0;
					for (int k = 0; k < j; ++k) {
						temp2 += pow(L.get(j, k), 2);
					}
					L.set(j, j, sqrt(temp - temp2));
				}

				if (i > j) {
					real temp = (1 / L.get(j, j));
					real temp2 = get(i, j);
					real temp3 = 0;
					for (int k = 0; k < j; ++k) {
						temp3 += (L.get(i, k) * L.get(j, k));
					}
					if (temp != 0) { L.set(i, j, temp * (temp2 - temp3)); }
				}

			}
		}
		return Matrix(rows, columns, L.element);
	}

	std::vector<Matrix> Matrix::LDL() {	//returns matrix L and D decomposition of matrix where A = LDL^T
		std::vector<Matrix> answer;
		if (rows <= 1 || columns <= 1) { return answer; }
		Matrix L(rows, columns);
		Matrix D(rows, columns);
		L.set(0, 0, sqrt(get(0, 0)));

		for (int j = 0; j < columns; ++j) {
			for (int i = j; i < rows; ++i) {
				if (i == j) {
					L.set(j, j, 1);
					real temp = 0;
					for (int k = 0; k < j; ++k) {
						temp += pow(L.get(j, k), 2) * D.get(k, k);
					}
					D.set(j, j, get(j, j) - temp);
				}

				if (i > j) {
					real temp = (1 / D.get(j, j));
					real temp2 = get(i, j);
					real temp3 = 0;
					for (int k = 0; k < j; ++k) {
						temp3 += (L.get(i, k) * L.get(j, k) * D.get(k, k));
					}
					L.set(i, j, temp * (temp2 - temp3));
				}
			}
		}

		answer.push_back(L);
		answer.push_back(D);
		return answer;
	}

	Matrix Matrix::dominantEigenvector() {//calculated by power method
		Matrix A(rows, columns, element);
		Matrix rd(A.rows, 1);
		rd.set(0, 0, 1);//set first value of matrix to 1
		Matrix init = rd;
		for (int i = 0; i < 500; ++i) {
			rd = A.multiply(A, rd);
			rd = rd.multiply(rd, 1 / rd.element[A.rows - 1]);
		}
		return rd;//dominant eigenvector
	}

	real Matrix::dominantEigenvalue() {//calculated by power method
		Matrix A = *this;
		Matrix rd(A.rows, 1);
		rd.set(0, 0, 1);//set first value of matrix to 1
		Matrix init = rd;
		for (int i = 0; i < 500; ++i) {
			rd = A.multiply(A, rd);
			rd = rd.multiply(rd, 1 / rd.element[A.rows - 1]);
		}
		real egv = A.multiply(A.multiply(A, rd), rd).sumColumn(0);
		egv /= rd.multiply(rd, rd).sumColumn(0);
		return egv;//dominant eigenvalue
	}

	bool Matrix::isPositiveDefinite() {//if all eigenvalues are positive, the matrix is positive definite
		Vector vec = eigenvaluesRealExact();
		for (int i = 0; i < vec.size(); ++i) {
			if (vec[i] <= 0) { return false; }
		}
		return true;
	}

	bool Matrix::isPositiveSemidefinite() {//if all eigenvalues are positive or 0, the matrix is positive semi-definite
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

	Vector Matrix::eigenvaluesRealExact() {//finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step
		Polynomial p = characteristicPolynomial();
		return p.realRoots();
	}

	std::vector<ComplexNumber> Matrix::eigenvaluesRealAndComplexExact() 
	{
		//Finds real eigenvalues by method of QR algorithm with Hessenberg intermediate step,
		//then either factors those roots to get the last few complex roots or uses the 
		//numerical Durand-Kerner method of root-finding

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
			if (std::abs(int(p.size()) - (int)realRoots.size()) <= 3) {
				//Use the shift method only if the real roots reduce it down to an easily calculable size (order 2 or less),
				//otherwise the results will be inaccurate
				if (p.size() <= 3) { i = realRoots.size(); break; }
				if (p.size() > 3) { p = p.factorOutBinomial(realRoots[i]); }
			}
		}
		if (p.size() < 3) { return complexRoots; }

		if (p.size() == 3) {
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
		{
			//start finding roots with Durand-Kerner method
			int iterations = 10000;
			ComplexNumber z = ComplexNumber(lowbound, lowboundC);
			int size = sizeof(z);
			std::vector<ComplexNumber> R;
			for (int i = 0; i < p.size(); i++) { 
				R.push_back(std::pow(z, i)); 
			}

			for (int i = 0; i < iterations; i++) {
				for (int j = 0; j < p.size(); j++) {
					ComplexNumber B = p.evaluate(R[j]);
					for (int k = 0; k < p.size(); k++) {
						if (k != j) { B /= (R[j] - R[k]); }
					}
					R[j] -= B;
				}
			}


			for (int i = 0; i < R.size(); ++i) { //now filter out all the bad results
				if (R[i].im() != 0) { //only complex roots accepted
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
							if (std::abs(a - x) < 0.001 && std::abs(b - y) < 0.001) { isSimilar = true; }
						}
						if (!isSimilar) { //if this is indeed a new root, save the root and its conjugate
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

	Vector Matrix::eigenvaluesNumerical() {//calculated by power method w/shifts
		Matrix A = *this;
		Matrix v;
		Vector egv;
		for (int i = 0; i < A.rows; ++i) {	//A' = A - (eigenvalue)*U*UT
			v = A.dominantEigenvector();
			real eigenval = A.multiply(A.multiply(A, v), v).sumColumn(0);
			real val = v.multiply(v, v).sumColumn(0);
			if (val != 0) { eigenval /= val; }
			else { eigenval = 0; }
			egv.push_back(eigenval);
			v.normalize();

			Matrix U = v.multiply(v, v.transpose());

			//U.display();
			U = U.multiply(U, -1 * eigenval);
			A = A.add(A, U);
		}
		return egv;
	}

	Matrix Matrix::eigenvalueMatrix() {
		Matrix A = *this;
		std::vector<Matrix> QR;

		int iterations = pow((A.rows * A.columns), 2);
		for (int i = 0; i < iterations; ++i) {//first, LU decomposition
			QR = A.LU();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}

		for (int i = 0; i < iterations; ++i) {//then, QR decomposition
			QR = A.QR();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}
		for (int i = 0; i < iterations; ++i) {//then, QR decomposition
			QR = A.HouseholderQR();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}

		return A;
	}

	Matrix Matrix::eigenvalueMatrix(Matrix A, int loops) {//recursive version of eigenvalueMatrix calculation
		if (loops <= 1) { return A; }
		--loops;
		std::vector<Matrix> QR;

		int iterations = (A.rows * A.columns);
		for (int i = 0; i < iterations; ++i) {//first, LU decomposition
			QR = A.LU();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}

		for (int i = 0; i < iterations; ++i) {//then, QR decomposition
			QR = A.QR();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}
		for (int i = 0; i < iterations; ++i) {//then, QR decomposition
			QR = A.HouseholderQR();
			if (QR[0].rows == 0 && QR[1].rows == 0) { break; }
			A = A.multiply(QR[1], QR[0]);
		}

		return eigenvalueMatrix(A, loops);
	}

	void Matrix::power(real n) {//exponentiates all elements in the matrix by some power
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				set(i, j, pow(get(i, j), n));
			}
		}
	}

	Matrix Matrix::leastSquares(Matrix X, Matrix Y) {//B=(X^T*X)^-1 * X^T * Y  
		Matrix XT = X.transpose();
		return Y.multiply(X.multiply(X.multiply(XT, X).inverse(), XT), Y);
	}

	Matrix Matrix::leastSquaresMatrix() {	//B=(X^T*X)^-1 * X^T * Y 
											/*note:  the matrix must be formatted in MINITAB format, such that the first column is the Y values, the rest are the
											column vectors of the X_i variable values. */

		Matrix M(rows, columns, element);//copy parent matrix
		if (columns != 2 && rows == 2) { M = M.transpose(); }//put matrix into column vector form

		Matrix X(M.rows, M.columns, element);	//isolate X variable values
		for (int i = 0; i < X.rows; ++i) {
			X.set(i, 0, 1);
		}

		Matrix Y(M.rows, 1, M.column(0));	//create seperate vector of y values

		return multiply(multiply(X.transpose(), X).inverse(), multiply(X.transpose(), Y));
	}

	Function Matrix::leastSquares() {	
		//B=(X^T*X)^-1 * X^T * Y 
		/*note:  the matrix must be formatted in MINITAB format, such that the first column is the Y values, the rest are the
		column vectors of the X_i variable values. See https://onlinecourses.science.psu.edu/stat501/node/292 for
		examples to prove this technique. There is test code also available at:  https://www.wessa.net/rwasp_multipleregression.wasp */

		Matrix M(rows, columns, element);//copy parent matrix
		if (columns != 2 && rows == 2) { M = M.transpose(); }//put matrix into column vector form

		Matrix X(M.rows, M.columns, element);//isolate X variable values
		for (int i = 0; i < X.rows; ++i) {
			X.set(i, 0, 1);
		}

		Matrix Y(M.rows, 1, M.column(0));//create seperate vector of y values

		Matrix coef = multiply(multiply(X.transpose(), X).inverse(), multiply(X.transpose(), Y));
		Vector cf;
		std::string funct = "f(";
		for (int i = 0; i < coef.rows; ++i) {
			cf.push_back(coef.get(i, 0));
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

	Vector Matrix::residuals() {
		//get least squares regression function
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

		//calculate resid. = y-yhat
		Vector residuals;
		for (int i = 0; i < yhat.size(); ++i) {
			residuals.push_back(get(i, 0) - yhat[i]);//in MINITAB format, column 0 is the y-value column
		}

		return residuals;
	}

	real Matrix::Rsquared() 
	{// aka Pearsson's "R"^2
	 //for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

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
			explainedError.push_back(pow(yhat[i] - ymean, 2));		// sum(yhat - yman)^2
			unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
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

	real Matrix::adjustedRsquared() {// aka Pearsson's "R"^2
									   //for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

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
			explainedError.push_back(pow(yhat[i] - ymean, 2));			// sum(yhat - yman)^2
			unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
		}

		//degress of freedom calculations
		//===============================
		real dfRegression = columns - 1;	//#predictor varibles aka 'k'
		real dfTotal = rows - 1;//total degrees freedom
		real dfResiduals = dfTotal - dfRegression;//# of residuals

		//sum of squares calculations
		real SSregression = explainedError.sum();
		real SSresiduals = unexplainedError.sum();
		real SStotal = totalError.sum();

		//mean square (MS) calculations
		//=============================
		real MSregression = SSregression / dfRegression;
		real MSresiduals = SSresiduals / dfResiduals;
		real MStotal = SStotal / dfTotal;

		//important statistics
		//====================
		real sqrtErrorVariance = sqrt(MSresiduals);	//equal to standard deviation of distribution of estimate
		real FTestStatistic = MSregression / MSresiduals;
		real Rsquared = (SStotal - SSresiduals) / SStotal;
		real adjustedRsquared = (MStotal - MSresiduals) / MStotal;

		return adjustedRsquared;
	}

	real Matrix::FTestStatistic() {
		//for explaination, see https://people.richland.edu/james/ictcm/2004/multiple.html

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
			explainedError.push_back(pow(yhat[i] - ymean, 2));			// sum(yhat - yman)^2
			unexplainedError.push_back(pow(get(i, 0) - yhat[i], 2));	// sum(y - yhat)^2
			totalError.push_back(pow(get(i, 0) - ymean, 2));			// sum(y - ymean)^2
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

	std::vector<Matrix> Matrix::thinQR() {
		std::vector<Matrix> answer = QR();
		answer[0];
		answer[1].trim();
		return answer;
	}

	std::vector<Matrix> Matrix::HouseholderQR() {//given some matrix, this produces its Householder QR factorization
		Matrix A = *this;
		Matrix Q(rows, columns);
		Q.identity();
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
			Matrix I(rows - j, columns - j);
			I.identity();
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

	Matrix Matrix::pseudoinverse() {//calculates the Moore-Penrose inverse: A* = R^-1 * Q^T.  Note:  A*A = I.
		std::vector<Matrix> qr = thinQR();
		return multiply(qr[1].inverse(), qr[0].transpose());
	}

	Polynomial Matrix::characteristicPolynomial() {
		Vector ans = characteristicMatrix().row(0);
		Vector answer(columns + 1);
		for (int i = 0; i < columns; ++i) { answer[columns - 1 - i] = -ans[i]; }
		answer[columns] = 1;
		return Polynomial(answer);
	}

	std::vector<Vector> Matrix::toVectorArray() {
		std::vector<Vector> arr;
		for (int i = 0; i < rows; ++i) {
			Vector ans;
			for (int j = 0; j < columns; ++j) {
				ans.push_back(get(i, j));
			}
			arr.push_back(ans);
		}
		return arr;
	}

	Matrix directionMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		A.identity();
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v[i]); }
		return A;
	}

	Matrix identityMatrix(int n) 
	{ 
		return identityMatrix(n, n); 
	}

	Matrix identityMatrix(int n, int m) {
		Vector vals(n * m);
		for (int i = 0; i < n; ++i) {
			for (int j = 0; j < m; ++j) {
				vals[i * m + j] = 0;
				if (i == j) { vals[i * m + j] = 1; }
			}
		}
		return Matrix(n, m, vals);
	}

	Matrix positionMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v[i]); }
		return A;
	}

	Matrix rotationMatrix(Vector v) {
		//Input rotation vector should be of form <angle,x,y,z> for 3D rotation, or <angle,x,y> for 2D.
		//If x, y, or z == 1, then it is a rotation around that axis.  Only angle and one coordinate should be non-zero
		if (v.size() < 3) { return Matrix(); }
		if (v.size() == 3) {//2D rotation matrix
			Matrix M(2, 2);
			M.set(0, 0, cos(v[0]));
			M.set(0, 1, -sin(v[0]));
			M.set(1, 0, sin(v[0]));
			M.set(1, 1, cos(v[0]));
			return M;
		}
		if (v.size() == 4) {//3D rotation matrix (contained in a 4x4 homogeneous coordinate matrix)
			if (v[1] == 1) {//x-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(1, 1, cos(v[0]));
				M.set(1, 2, -sin(v[0]));
				M.set(2, 1, sin(v[0]));
				M.set(2, 2, cos(v[0]));
				return M;
			}
			if (v[2] == 1) {//y-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(0, 0, cos(v[0]));
				M.set(0, 2, -sin(v[0]));
				M.set(2, 0, sin(v[0]));
				M.set(2, 2, cos(v[0]));
				return M;
			}
			if (v[3] == 1) {//z-axis rotation
				Matrix M(4, 4);
				M.identity();
				M.set(0, 0, cos(v[0]));
				M.set(0, 1, -sin(v[0]));
				M.set(1, 0, sin(v[0]));
				M.set(1, 1, cos(v[0]));
				return M;
			}
		}
		return Matrix();
	}

	Matrix scalingMatrix(Vector v) {
		int sz = v.size() + 1;
		Matrix A(sz, sz);
		A.identity();
		for (int i = 0; i < v.size(); ++i) { A.set(i, i, v[i]); }
		return A;
	}

	Matrix translationMatrix(Vector v) {
		if (v.size() <= 1) { return Matrix(); }
		if (v.size() == 2) {}
		if (v.size() == 3) {
			Matrix A(4, 4);
			A.identity();
			for (int i = 0; i < v.size(); ++i) { A.set(i, A.columns - 1, v[i]); }
			return A;
		}
		return Matrix();
	}

	Matrix companionMatrix(Polynomial p) {//Froebenius form
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

	Matrix Vandermonde(Polynomial p) {
		//create a Vandermonde matrix from the vector of coefficients from a_n-1 to a0,
		//which is suitable for finding the complex eigenvalues of the polynomial
		Matrix CM = companionMatrix(p);
		Vector co(p);
		Matrix V((p.size() - 1), p.size() - 1);
		for (int i = 0; i < p.size() - 1; ++i) {
			for (int j = 0; j < p.size() - 1; ++j) {
				if (i == 0) { V.set(i, j, 1); }
				else { V.set(i, j, pow(co[i], j + 1)); }
			}
		}
		return V;
	}
}