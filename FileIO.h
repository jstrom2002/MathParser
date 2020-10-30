/**
*	Header for file I/O, mostly for use with loading grayscale images of various types to
*	'Matrix' objects.
*/

#pragma once
#include "Types.h"
#include <vector>

namespace MathParser
{
	class Matrix;

	Matrix csvread(std::string filename);
	void csvwrite(std::string filename, Matrix& M);

	// Functions for loading color image data from file to an array of grayscale reals.
	void loadBMP(std::string filename, int* rows, int* cols, std::vector<real>* pixels);
	void loadPPM(std::string filename, int* rows, int* cols, std::vector<real>* pixels);
	void loadTGA(std::string filename, int* rows, int* cols, std::vector<real>* pixels);
	void saveBMP(std::string filename, int rows, int cols, int colorDepth, 
		std::vector<real>* pixels);
	void savePPM(std::string filename, int rows, int cols, int colorDepth,
		std::vector<real>* pixels);
	void saveTGA(std::string filename, int rows, int cols, int colorDepth, 
		std::vector<real>* pixels);
}