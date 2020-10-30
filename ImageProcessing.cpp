#include "ImageProcessing.h"
#include "Matrix.h"
#include "Vector.h"
#include "StringUtils.h"
#include <fstream>
#include "Vector.h"

namespace MathParser
{
	real convertRGBtoGray(Vector rgb)
	{	// This function is intended to handle byte-valued data from loaded objects,
		// hence normalization by 255. 'coef' values are gamma-corrected luminance 
		// conversion values.
		const real coef[3] = { 0.299 / 255.0, 0.587 / 255.0, 0.114 / 255.0 };
		return (coef[0]*rgb[0] + coef[1]*rgb[1] + coef[2]*rgb[2]);
	}

	Vector convertGrayToRGB(real gray)
	{
		// For now, this function simply returns a 3-vector of the grayscale value
		// clamped to valid byte-range [0,255].
		real R = gray;
		real G = gray;	
		real B = gray;
		Vector rgb(std::vector<real>{ R, G, B });
		rgb *= 255;
		rgb = clamp(rgb, 0, 255);
		return rgb;
	}

	Matrix identityFilter(unsigned int sz)
	{	// Filter which, when convolved with a matrix, should not alter the image in any way.
		if (sz % 2 == 0) //filters must be of odd size
			return Matrix();

		std::vector<real> vals(sz * sz, 0);
		vals[std::floor((sz * sz) / 2)] = 1;
		Matrix Mat(sz, sz, vals);
		return Mat;
	}

	Matrix sharpeningFilter(unsigned int sz) {
		if (sz % 2 == 0) //filters must be of odd size
			return Matrix();

		std::vector<real> vals;
		if (sz == 3) {
			vals = std::vector<real>{
				 0, -1,  0,
				-1,  5, -1,
				 0, -1,  0
			};
		}
		Matrix Mat(sz, sz, vals);
		return Mat;
	}

	Matrix edgeDetectorFilter(unsigned int sz) {
		if (sz % 2 == 0) //filters must be of odd size
			return Matrix();

		std::vector<real> vals;
		if (sz == 3) {
			vals = std::vector<real>{
				-1, -1, -1,
				-1,  8, -1,
				-1, -1, -1
			};
		}

		Matrix Mat(sz, sz, vals);
		return Mat;
	}

	Matrix boxFilter(unsigned int sz)
	{
		if (sz % 2 == 0) { return Matrix(); }//filters must be of odd size
		std::vector<real> vals(sz * sz, 1);
		vals[std::floor((sz * sz) / 2) + 1] = 1;
		Matrix Mat(sz, sz, vals);
		return Mat;
	}

	Matrix GaussianKernel1D(int sz, real sigma)
	{
		real sigmaX = sigma > 0 ? sigma : ((sz - 1) * 0.5 - 1) * 0.3 + 0.8;
		real scale2X = -0.5 / (sigmaX * sigmaX);
		real sum = 0;

		std::vector<real> cd(sz, 0);

		// For smaller sized kernels, use precalculated array of values.
		const int SMALL_GAUSSIAN_SIZE = 7;
		static const float small_gaussian_tab[][SMALL_GAUSSIAN_SIZE] = {
			{1.f},
			{0.25f, 0.5f, 0.25f},
			{0.0625f, 0.25f, 0.375f, 0.25f, 0.0625f},
			{0.03125f, 0.109375f, 0.21875f, 0.28125f, 0.21875f, 0.109375f, 0.03125f}
		};

		// For larger sized kernels, do calculation.
		const float* fixed_kernel = sz % 2 == 1 && sz <= SMALL_GAUSSIAN_SIZE &&
			sigma <= 0 ? small_gaussian_tab[sz >> 1] : 0;
		int i;
		for (i = 0; i < sz; i++) {
			real x = i - (sz - 1) * 0.5;
			real t = fixed_kernel ? (real)fixed_kernel[i] : std::exp(scale2X * x * x);
			cd[i] = t;
			sum += cd[i];
		}
		sum = 1.0 / sum;
		for (i = 0; i < sz; i++) { cd[i] *= sum; }
		return Matrix(sz, 1, cd);
	}

	Matrix GaussianKernel2D(int sz, real sigma)
	{
		Matrix Mat1 = GaussianKernel1D(sz, sigma);
		Matrix MatT = Mat1.transpose();
		return Mat1 * MatT;
	}

	Matrix gradientFilter(std::string type, bool dy)
	{
		int sz1, sz2;
		std::vector<real> vals;
		if (type == "finite difference" || type == "") {
			vals = std::vector<real>{ 0,-1,1 };
			sz1 = 3;
			sz2 = 1;
		}
		if (type == "Sobel") {
			std::vector<real> vals = std::vector<real>{
				 -1, 0, 1,
				 -2, 0, 2,
				 -1, 0, 1
			};
			sz1 = 3;
			sz2 = 3;
		}
		if (type == "Schurr") {
			std::vector<real> vals = std::vector<real>{
				 -3, 0, 3,
				-10, 0,10,
				 -3, 0, 3
			};
			sz1 = 3;
			sz2 = 3;
		}
		if (type == "Laplacian") {
			std::vector<real> vals = std::vector<real>{
				0, 1, 0,
				1,-4, 1,
				0, 1, 0
			};
			sz1 = 3;
			sz2 = 3;
		}
		Matrix Mat(sz1, sz2, Vector(vals));
		if (dy)
			Mat = Mat.transpose();
		return Mat;
	}
}