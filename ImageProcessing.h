#pragma once
#include "Types.h"
#include "BaseObject.h"
#include <string>
#include <vector>

namespace MathParser
{
	class Matrix;
	class Vector;

	class IMAGE : public BaseObject
	{
	public:
		int colorDepth; //read color depth to determine color array
		int compressionMethod;//tells what method of compression is used (0 = uncompressed) 
		real gamma;//gamma value for the immage, Vout = A*Vin^(gamma)
		int headerSize;//size of header before pixel array
		int hgt;//gives image width
		int horizontalResolution;
		int imageSize;//gives image size
		int numberOfBytes;//read number of bytes in file
		int numberOfColorsInPalatte;//colors in palatte, 0 is default (2^n colors)
		int numberOfPixels;//number of total pixel vectors, ie for 24-bit arrays it would be 3 * wdt * hgt
		int pixelsPerRow;//number of pixels in a row, width*3
		int padBytes;//number of bytes padding each row
		int reservedBytes;//field of the header
		int rowSize;//number of bytes necessary to store one row
		int verticalResolution;
		int wdt;//gives image height 
		std::string filename;
		std::string filetype;//the extension of the image file (ie BMP, JPEG, GIF, etc)
		std::string headertype;
		std::string colorOrdering;//displays whether the array values are RGB, BGR, grayscale, etc
		std::vector<unsigned char> header;//byte vector of the data for the header to BMP file
		std::vector<std::vector<unsigned char>> pixelArray;//an array of all the pixels in the BMP file 

		IMAGE() {}
		IMAGE(std::string fn);
		IMAGE(Matrix A, std::string fileType="BMP");
		void convertToGrayscale();
		Matrix getBlackAndWhiteMatrix();
		void saveImage(std::string name);
		void printToTextFile(std::string outputFile);
		Vector toVector();
		std::vector<unsigned char> get(int i, int j);
		unsigned char get(int i, int j, int n);
		void set(int i, int j, std::vector<unsigned char> val);
		void set(int i, int j, int n, unsigned char val);
		Matrix getRedMatrix();
		Matrix getGreeMatrix();
		Matrix getBlueMatrix();
		void maskOutRed();
		void maskOutGreen();
		void maskOutBlue();
		void invertColor();
		void transpose();
		void convertBGRtoRBG();
		void adjustGamma(real gam);
		std::vector<unsigned char> makeHeader(std::string fileType);
		void updateHeader();
		int calculateRowSize(std::string type);
		int calculateNumberOfBytes(std::string type);
		std::string getHeadertype();
		virtual std::string to_string(int precision = 4);

	private:
		unsigned long* make_crc_table();

	};
	bool isBMP(std::string filename);
	bool isGIF(std::string filename);
	bool isJPEG(std::string filename);
	bool isPNG(std::string filename);

	void fileToText(std::string filename);
	void BMPtoText(std::string filename, std::string outputName="HexDump.txt");

	Matrix identityFilter(unsigned int sz);
	Matrix sharpeningFilter(unsigned int sz);
	Matrix edgeDetectorFilter(unsigned int sz);
	Matrix boxFilter(unsigned int sz);
	Matrix GaussianKernel1D(int sz, real sigma=1);
	Matrix GaussianKernel2D(int sz, real sigma=1);
	Matrix gradientFilter(std::string type, unsigned int wrt=0);
	Matrix testMatrix();
}