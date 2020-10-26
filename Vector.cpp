#include "Vector.h"
#include "Matrix.h"
#include "Polynomial.h"
#include <limits>
#include <iomanip>

namespace std
{
	MathParser::Vector pow(MathParser::Vector v, int n)
	{
		if (!n)//Handle case where exponent == 0.
		{
			std::vector<MathParser::real> vals;
			vals.resize(v.size(),1);
			return MathParser::Vector(vals);
		}

		MathParser::Vector v_2 = v;
		for (int i = 0; i < std::abs(n); ++i)
			v_2 *= v;

		if (n < 0)//Handle negative exponents.
			return 1.0 / v_2;
		else
			return v_2;
	}
}

namespace MathParser
{
	Vector::Vector(std::vector<real> v) : v(v)
	{
	}

	Vector::Vector(Polynomial p) 
	{
		resize(p.size(),0);
		for (int i = 0; i < p.size(); ++i)
			v[i] = p[i];
	}

	void Vector::removeNullValues()
	{
		for (int i = 0; i < v.size(); ++i)
			if (std::abs(v[i]) < 0) 
				v.erase(v.begin() + i); 
	}

	Matrix Vector::toDiagonalMatrix() 
	{
		Matrix M(v.size(), v.size());
		for (int i = 0; i < size(); ++i) 		
			M.set(i, i, v[i]);		
		return M;
	}

	Matrix Vector::toMatrix() 
	{ 
		return Matrix(1, size(), v); 
	}

	Matrix Vector::toMatrix(int r, int c) 
	{
		Matrix M(r, c);
		M.identity();
		for (int i = 0; i < size(); ++i) 
			M.set(0, i, v[i]);		
		return M;
	}

	std::string Vector::to_string(int precision) 
	{
		int vectorPrecision = precision;
		std::ostringstream sstr;
		sstr << "(";
		for (int i = 0; i < v.size(); ++i) {
			if (v[i] < std::numeric_limits<real>::max() && 
				v[i] > -std::numeric_limits<real>::max()) {
				if (v[i] == floor(v[i])) 
					vectorPrecision = 0;
				else 
					vectorPrecision = precision;
				sstr << std::setprecision(vectorPrecision) << v[i];
				if (i < v.size() - 1) { sstr << ","; }
			}
		}
		sstr << ")";
		std::string answer = sstr.str();
		sstr.clear();
		return answer;
	}

	real dot(Vector v1, Vector v2)
	{
		real val = 0;
		unsigned int sz = v1.size() < v2.size() ? v1.size() : v2.size();
		for (unsigned int i = 0; i < sz; ++i)
		{
			val += v1[i] * v2[i];
		}
		return val;
	}

	Vector cross(Vector a, Vector b)
	{
		int sz = a.size() < b.size() ? a.size() : b.size();
		// Depending upon number of dimensions of vector, 
		// cross product is calculable by either simple formula
		// or by more complex Hodge star operator/wedge product.
		if (sz != 3)
			return Vector();

		return Vector(std::vector<real>{
			a[1]*b[2]-a[2]*b[1],
			real(-1.0)*(a[0]*b[2]-a[2]*b[0]),
			a[0]*b[1]-a[1]*b[0]
		});
	}

	real Vector::length()
	{
		real sum = 0;
		for (int i = 0; i < v.size(); ++i)
			sum += v[i] * v[i];
		sum = std::sqrt(sum);
		return sum;
	}

	real length(Vector v)
	{
		real sum = 0;
		for (int i = 0; i < v.size(); ++i)
			sum += v[i] * v[i];
		sum = sqrt(sum);
		return sum;
	}

	Vector normalize(Vector v)
	{
		real n = v.length();
		for (int i = 0; i < v.size(); ++i)
			v[i] /= n;
		return v;
	}

	real Vector::mean()
	{
		return sum() / size();
	}

	real Vector::sum()
	{
		real m = 0;
		for (int i = 0; i < v.size(); ++i)
			m += v[i];
		return m;
	}

	real Vector::sumexp(real p)
	{
		real m = 0;
		for (int i = 0; i < v.size(); ++i)
			m += pow(v[i],p);
		return m;
	}

	real Vector::populationVariance() 
	{
		real mu = mean();
		real answer = 0;
		for (int i = 0; i < size(); ++i) {
			answer += pow(v[i] - mu, 2);
		}
		return answer / size();
	}

	real Vector::populationStandardDeviation() 
	{
		return sqrt(populationVariance());
	}

	real Vector::sampleVariance() 
	{
		int sz = size();
		real mu = mean();
		real answer = 0;
		for (int i = 0; i < sz; ++i) 
		{
			answer += pow(v[i] - mu, 2);
		}
		return answer / (sz - 1);
	}

	real Vector::sampleStandardDeviation() 
	{
		return sqrt(sampleVariance());
	}

	real Vector::RMS() 
	{
		real answer = 0;
		for (int i = 0; i < size(); ++i) {
			answer += pow(v[i], 2);
		}
		answer *= 1.0 / size();
		return pow(answer, 0.5);
	}

	Vector Vector::add(Vector a, Vector b)
	{
		std::vector<real> v;
		int small_sz = a.size() < b.size() ? a.size() : b.size();
		int large_sz = a.size() > b.size() ? a.size() : b.size();
		v.resize(large_sz);
		for (int i = 0; i < large_sz; ++i)
		{
			if (i < small_sz)
				v[i] = a[i] + b[i];
			else if (a.size() == large_sz)
				v[i] = a[i];
			else if (b.size() == large_sz)
				v[i] = b[i];
		}
		return Vector(v);
	}

	Vector Vector::add(Vector a, real  b)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{			
			v[i] = a[i] + b;
		}
		return Vector(v);
	}
	
	Vector Vector::add(real b, Vector a)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] + b;
		}
		return Vector(v);
	}

	Vector Vector::subtract(Vector a, Vector b)
	{
		std::vector<real> v;
		int small_sz = a.size() < b.size() ? a.size() : b.size();
		int large_sz = a.size() > b.size() ? a.size() : b.size();
		v.resize(large_sz);
		for (int i = 0; i < large_sz; ++i)
		{
			if (i < small_sz)
				v[i] = a[i] - b[i];
			else if (a.size() == large_sz)
				v[i] = a[i];
			else if (b.size() == large_sz)
				v[i] = b[i];
		}
		return Vector(v);
	}

	Vector Vector::subtract(Vector a, real  b)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] - b;
		}
		return Vector(v);
	}

	Vector Vector::subtract(real b, Vector a)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] - b;
		}
		return Vector(v);
	}

	Vector Vector::multiply(Vector a, Vector b)
	{
		std::vector<real> v;
		int small_sz = a.size() < b.size() ? a.size() : b.size();
		int large_sz = a.size() > b.size() ? a.size() : b.size();
		v.resize(large_sz);
		for (int i = 0; i < large_sz; ++i)
		{
			if (i < small_sz)
				v[i] = a[i] * b[i];
			else if (a.size() == large_sz)
				v[i] = a[i];
			else if (b.size() == large_sz)
				v[i] = b[i];
		}
		return Vector(v);
	}

	Vector Vector::multiply(Vector a, real b)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] * b;
		}
		return Vector(v);
	}

	Vector Vector::multiply(real b, Vector a)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] * b;
		}
		return Vector(v);
	}

	Vector Vector::divide(Vector a, Vector b)
	{
		std::vector<real> v;
		int small_sz = a.size() < b.size() ? a.size() : b.size();
		int large_sz = a.size() > b.size() ? a.size() : b.size();
		v.resize(large_sz);
		for (int i = 0; i < large_sz; ++i)
		{
			if (i < small_sz)
				v[i] = a[i] / b[i];
			else if (a.size() == large_sz)
				v[i] = a[i];
			else if (b.size() == large_sz)
				v[i] = b[i];
		}
		return Vector(v);
	}

	Vector Vector::divide(Vector a, real b)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] / b;
		}
		return Vector(v);
	}
	
	Vector Vector::divide(real b, Vector a)
	{
		std::vector<real> v;
		int sz = a.size();
		v.resize(sz);
		for (int i = 0; i < sz; ++i)
		{
			v[i] = a[i] / b;
		}
		return Vector(v);
	}

	Vector& Vector::operator*=(Vector& rhs) { *this = *this * rhs; return *this; }
	Vector& Vector::operator*=(real x) { *this = *this * x; return *this; }
	Vector& Vector::operator/=(Vector& rhs) { *this = *this / rhs; return *this; }
	Vector& Vector::operator/=(real x) { *this = *this / x; return *this; }
	Vector& Vector::operator+=(Vector& rhs) { *this = *this + rhs; return *this; }
	Vector& Vector::operator+=(real x) { *this = *this + x; return *this; }
	Vector& Vector::operator-=(Vector& rhs) { *this = *this - rhs; return *this; }
	Vector& Vector::operator-=(real x) { *this = *this - x; return *this; }
	Vector operator+(Vector lhs, Vector rhs) { return lhs.add(lhs, rhs); }
	Vector operator+(Vector lhs, real x) { return lhs.add(lhs, x); }
	Vector operator+(real x, Vector rhs) { return rhs.add(x, rhs); }
	Vector operator-(Vector lhs, Vector rhs) { return lhs.subtract(lhs, rhs); }
	Vector operator-(Vector lhs, real x) { return lhs.subtract(lhs, x); }
	Vector operator-(real x, Vector rhs) { return rhs.subtract(x, rhs); }
	Vector operator*(Vector lhs, Vector rhs) { return lhs.multiply(lhs, rhs); }
	Vector operator*(Vector lhs, real x) { return lhs.multiply(lhs, x); }
	Vector operator*(real x, Vector rhs) { return rhs.multiply(x, rhs); }
	Vector operator/(Vector lhs, Vector rhs) { return lhs.divide(lhs, rhs); }
	Vector operator/(Vector lhs, real x) { return lhs.divide(lhs, x); }
	Vector operator/(real x, Vector rhs) { return rhs.divide(x, rhs); }
}