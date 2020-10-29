#include "FileIO.h"
#include "ImageProcessing.h"
#include "Parsing.h"
#include "Matrix.h"
#include "StringUtils.h"
#include <fstream>
#include <iostream>
#include <string>
#include <cmath>
#include <algorithm>

namespace MathParser
{
	Matrix csvread(std::string filename)
	{
		std::vector<real> element;
		std::ifstream file1(filename, std::ios::in);
		std::string buff = "";
		int r = 0;
		int c = 0;
		while (std::getline(file1, buff))
		{
			while (buff[0] < '0' || buff[0] > 'z')
				buff = buff.substr(1);
			while (buff[buff.length() - 1] < '0' || buff[buff.length() - 1] > 'z')
				buff = buff.substr(0, buff.length() - 1);
			std::vector<real> temp = ParseNInputNumbers(buff);
			if (temp.size() > c)
				c = temp.size();
			element.insert(std::end(element), std::begin(temp), std::end(temp));
			r++;
		}
		return Matrix(r, c, element);
	}

	void csvwrite(std::string filename, Matrix& M)
	{
		FILE* file1 = fopen(filename.c_str(), "w");
		while (file1)
		{
			for (int i = 0; i < M.rows; ++i)
			{
				for (int j = 0; j < M.columns; ++j)
				{
					if (j < M.columns - 1)
						printf("%d,", M.element[i * M.columns + j]);
					else if (i < M.rows - 1)
						printf("%d\n", M.element[i * M.columns + j]);
					else
						printf("%d", M.element[i * M.columns + j]);
				}
			}
			fclose(file1);
		}
	}

	void loadBMP(std::string filename, int* rows, int* cols, std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "rb");
		if (!f)//If file does not exist, return.
			return;

		int numberOfBytes = 0;
		int reservedBytes = 0;
		int headerSize = 0;
		int colorDepth = 0;

		//read preliminary file data -- 14 bytes
		std::vector<unsigned char> header;
		unsigned char prelimData[14];
		fread(prelimData, sizeof(unsigned char), 14, f);
		for (int i = 0; i < 14; ++i)
			header.push_back(prelimData[i]);
		numberOfBytes = (header[5] << 24) ^ (header[4] << 16) ^ (header[3] << 8) ^ header[2];//read number of bytes in file
		reservedBytes = (header[9] << 24) ^ (header[8] << 16) ^ (header[7] << 8) ^ header[6];//read reserved data
		headerSize = (header[13] << 24) ^ (header[12] << 16) ^ (header[11] << 8) ^ header[10];//read starting address

		//read and interpret file header data
		unsigned char* headerData = new unsigned char[headerSize - 14];
		fread(headerData, sizeof(unsigned char), headerSize - 14, f);//read the 54-byte header
		for (int i = 0; i < headerSize - 14; ++i)
			header.push_back(headerData[i]);
		delete[] headerData;

		//initialize class variables;
		int DIBheaderSize = (header[17] << 24) ^ (header[16] << 16) ^ (header[15] << 8) ^ header[14];
		*cols = (header[21] << 24) ^ (header[20] << 16) ^ (header[19] << 8) ^ header[18];
		*rows = (header[25] << 24) ^ (header[24] << 16) ^ (header[23] << 8) ^ header[22];
		int numberOfPlanes = (header[27] << 8) ^ header[26];
		colorDepth = (header[29] << 8) ^ header[28]; //read color depth to determine color array
		int compression = (header[33] << 24) ^ (header[32] << 16) ^ (header[31] << 8) ^ header[30];
		int dataSizeWPadding = (header[37] << 24) ^ (header[36] << 16) ^ (header[35] << 8) ^ header[34];
		int pixelsPerMeterX = (header[41] << 24) ^ (header[40] << 16) ^ (header[39] << 8) ^ header[38];
		int pixelsPerMeterY = (header[45] << 24) ^ (header[44] << 16) ^ (header[43] << 8) ^ header[42];
		int numberColorsInPalatte = (header[49] << 24) ^ (header[48] << 16) ^ (header[47] << 8) ^ header[46];
		int significantColors = (header[53] << 24) ^ (header[52] << 16) ^ (header[51] << 8) ^ header[50];

		// Find out # of padding bits in a row.
		int padBits = (4 - ((*cols * 3) % 4)) % 4;
		int rowSize = *cols + padBits / (colorDepth / 8);

		// Now that the header has been read, load all subsequent bytes.
		std::vector<unsigned char> data(dataSizeWPadding, 0);
		fread(data.data(), sizeof(unsigned char), dataSizeWPadding, f);
		fclose(f);

		if (colorDepth == 8)// If already grayscale, just copy array of data to element array.
		{
			for (int i = 0; i < data.size(); ++i)
				pixels->push_back(data[i]);
		}
		else// Else, convert to grayscale.
		{
			int pixelCount = 0;//counts the number of pixel vectors read
			std::vector<unsigned char> tempVec;
			if (colorDepth == 24)
			{
				for (int i = 0; i < data.size(); ++i) {
					if ((i % rowSize) < (rowSize - padBits / (colorDepth / 8)))
					{
						tempVec.push_back(data[i]);
					}
					if (tempVec.size() >= 3 && i % 3 == 0)
					{	// Calculate luminance to convert to grayscale.
						std::vector<real> vals = std::vector<real>{ real(tempVec[2]),
							real(tempVec[1]), real(tempVec[0]) };
						pixels->push_back(convertRGBtoGray(vals));
						tempVec.clear();
					}
				}
			}
		}
	}

	void loadPPM(std::string filename, int* rows, int* cols, std::vector<real>* pixels)
	{
		std::ifstream filein(filename);
		std::string buff = "";
		int maxValue = 0;
		int headerSize = 0;
		std::vector<std::string> strs;

		// Get all header info first using ifstream.
		getline(filein,buff);
		headerSize += buff.length()+1;
		if (buff.find("P6") == std::string::npos && buff.find("P3") == std::string::npos)
		{
			filein.close();
			return;
		}
		// Next values are columns and rows.
		getline(filein, buff);
		headerSize += buff.length()+1;
		while (buff[0] == '#' || !buff.length())//ignore comments in file.
		{
			getline(filein, buff);
			headerSize += buff.length()+1;
		}
		strs = tokenize(buff, " ");
		if(strs.size() > 0)
			*cols = std::stoi(strs[0]);
		if (strs.size() > 1)
			*rows = std::stoi(strs[1]);

		// Next value is maximum pixel value.
		getline(filein, buff);
		headerSize += buff.length()+1;
		while (buff[0] == '#' || !buff.length())//ignore comments in file.
		{
			getline(filein, buff);
			headerSize += buff.length()+1;
		}
		maxValue = std::stoi(buff);
		filein.close();

		// Now read pixels using a FILE pointer.
		FILE* f = fopen(filename.c_str(), "rb");
		if (!f)//If file does not exist, return.
			return;
		fseek(f, headerSize, SEEK_CUR);//skip past header.

		int sz = 3 * (*cols) * (*rows);
		std::vector<unsigned char> px(sz, 0);
		fread(px.data(), 1, sz, f);
		for(int i=0; i<sz; ++i)
		{
			std::vector<real> rls = std::vector<real>{ real(px[i+0]),
				real(px[i+1]), real(px[i+2]) };
			pixels->push_back(convertRGBtoGray(rls));
		}
		fclose(f);
	}

	void loadTGA(std::string filename, int* rows, int* cols, std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "rb");
		if (!f)//If file does not exist, return.
			return;

		// Get info from header.
		unsigned char TGAheader[12];
		unsigned char header[6];
		fread(TGAheader, sizeof(unsigned char), 12, f);
		fread(header, sizeof(unsigned char), 6, f);
		*cols = (header[1] << 8) ^ header[0];
		*rows = (header[3] << 8) ^ header[2];
		int colorDepth = header[4];

		// Now read pixel data.
		std::vector<unsigned char> temp;
		int nSize = (*rows) * (*cols) * (colorDepth / 8);
		temp.resize(nSize, 0);
		fread(temp.data(), sizeof(unsigned char), nSize, f);

		// Copy pixel data to array.
		for (int i = 0; i < nSize - (colorDepth / 8); i += (colorDepth / 8))
		{
			if ((colorDepth / 8) == 1)
				pixels->push_back(temp[i]);
			else if ((colorDepth / 8) == 3 || (colorDepth / 8) == 4)
			{
				std::vector<real> vec = std::vector<real>{ real(temp[i + 2]), 
					real(temp[i + 1]), real(temp[i + 0]) };
				pixels->push_back(convertRGBtoGray(vec));
			}
		}
	}

	void saveBMP(std::string filename, int rows, int cols, int colorDepth,
		std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "wb");
		int padBits = (4 - ((cols * 3) % 4)) % 4;
		int padPixels = padBits / (colorDepth / 8);//assumes colorDepth is bits per pixel, ie 24-bit.

		// Write first 14 byte BMP ID header.
		uint32_t fileBytes = (cols + padPixels) * rows + 14 + 40;
		uint32_t offsetToData = 14 + 40;
		uint32_t dummy = 0;
		unsigned char BMPheader[14] = {
			0x42, 0x4D,
			(fileBytes), (fileBytes >> 8), (fileBytes >> 16), (fileBytes >> 24),
			0,0,0,0,
			(offsetToData), (offsetToData >> 8), (offsetToData >> 16), (offsetToData >> 24)
		};
		fwrite(BMPheader, sizeof(unsigned char), 14, f);

		// Write 40 byte DIB header.
		uint32_t bmpSize = rows * (cols + padPixels);
		unsigned char DIBheader[40] =
		{
			40, 0, 0, 0,
			(cols), (cols >> 8), (cols >> 16), (cols >> 24),
			(rows), (rows >> 8), (rows >> 16), (rows >> 24),
			1, 0, colorDepth, 0,
			0, 0, 0, 0,
			(bmpSize), (bmpSize >> 8), (bmpSize >> 16), (bmpSize >> 24),
			(2834), (2834 >> 8), (2834 >> 16), (2834 >> 24),
			(2834), (2834 >> 8), (2834 >> 16), (2834 >> 24),
			0, 0, 0, 0,
			0, 0, 0, 0
		};
		fwrite(DIBheader, sizeof(unsigned char), 40, f);

		//write pixel data.
		unsigned char val0 = 0;
		int counter = 0;
		for (int i = 0; i < rows; i++)
		{
			for (int j = 0; j <= cols; j++)
			{
				if (j < cols)
				{
					// Convert back from Gray to RGB and write pixels to file.
					real val = (*pixels)[i * cols + j];
					std::vector<real> rgb = convertGrayToRGB(val);
					fwrite(rgb.data(), sizeof(unsigned char), 3, f);
					counter += 3;
				}
				else
				{
					for (int n = 0; n < padBits; ++n)
						fwrite(&val0, 1, 1, f);
					counter += padBits;
				}
			}
		}

		// Write values to array as necessary.
		while (counter < (rows * (cols + padPixels)))
		{
			fwrite(&val0, 1, 1, f);
			counter++;
		}
		fclose(f);
	}

	void savePPM(std::string filename, int rows, int cols, int colorDepth,
		std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "wb");

		// Somewhat clumsy gathering of header values into a char array for writing.
		std::vector<char> PPMheader = { 'P', '6', '\n' };
		
		std::string colstr = std::to_string(cols) + " ";
		std::vector<char> colstr2(colstr.begin(), colstr.end());
		//std::reverse(colstr2.begin(), colstr2.end());//reverse byte order.

		std::string rowstr = std::to_string(rows) + "\n";
		std::vector<char> rowstr2(rowstr.begin(), rowstr.end());
		//std::reverse(rowstr2.begin(), rowstr2.end());//reverse order of bytes.

		std::vector<char> maxval{'2', '5', '5'};
		PPMheader.insert(std::end(PPMheader), std::begin(colstr2), std::end(colstr2));
		PPMheader.insert(std::end(PPMheader), std::begin(rowstr2), std::end(rowstr2));
		PPMheader.insert(std::end(PPMheader), std::begin(maxval), std::end(maxval));

		// Write header.
		fwrite(PPMheader.data(), sizeof(unsigned char), PPMheader.size(), f);

		// Convert reals to bytes. PPM files are only RGB valued pixels.
		size_t sz = rows * cols;
		std::vector<unsigned char> temp(sz * 3, 0);
		for (int i = 0; i < pixels->size(); ++i)
			temp[i] = (*pixels)[i];
		fwrite(temp.data(), sizeof(unsigned char), sz * 3, f);
		fclose(f);
	}

	void saveTGA(std::string filename, int rows, int cols, int colorDepth,
		std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "wb");
		int nSize = cols * rows * 3;

		// Transfer elements to pixel array.
		for (int i = 0; i < rows; i++)
			for (int j = 0; j < cols; j++)
			{
				// Convert back from Gray to RGB and write pixels to file.
				real val = (*pixels)[i * cols + j];
				std::vector<real> rgb = convertGrayToRGB(val);
				pixels->push_back(rgb[2]);//convert to BGR.
				pixels->push_back(rgb[1]);
				pixels->push_back(rgb[0]);
			}

		// Make headers, write all to file.
		unsigned char TGAheader[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
		unsigned char header[6] = { cols % 256, cols / 256,
			rows % 256, rows / 256, colorDepth, 0 };
		fwrite(TGAheader, sizeof(unsigned char), 12, f);
		fwrite(header, sizeof(unsigned char), 6, f);
		fwrite(pixels->data(), sizeof(unsigned char), nSize, f);	
		fclose(f);
	}
}