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
#include <array>

namespace MathParser
{
	void MergeBytes(std::array<unsigned char, 4>* pixel,
		unsigned char* p, int bytes)
	{
		if (bytes == 4) {
			(*pixel)[0] = p[2];
			(*pixel)[1] = p[1];
			(*pixel)[2] = p[0];
			(*pixel)[3] = p[3];
		}
		else if (bytes == 3) {
			(*pixel)[0] = p[2];
			(*pixel)[1] = p[1];
			(*pixel)[2] = p[0];
			(*pixel)[3] = 255;
		}
		else if (bytes == 2) {
			(*pixel)[0] = (p[1] & 0x7c) << 1;
			(*pixel)[1] = ((p[1] & 0x03) << 6) | ((p[0] & 0xe0) >> 2);
			(*pixel)[2] = (p[0] & 0x1f) << 3;
			(*pixel)[3] = (p[1] & 0x80);
		}
	}


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
		if (!file1)
			return;

		for (int i = 0; i < M.rows; ++i)
		{
			for (int j = 0; j < M.columns; ++j)
			{
				if (j < M.columns - 1)
					fprintf(file1, "%f,", M(i,j));
				else if (i < M.rows - 1)
					fprintf(file1, "%f\n", M(i, j));
				else
					fprintf(file1, "%f", M(i, j));
			}
		}
		
		fclose(file1);
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
						std::vector<real> vals = std::vector<real>{ real(tempVec[0]),
							real(tempVec[1]), real(tempVec[2]) };
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
		for(int i=0; i<sz; i+=3)
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
		
		/*  TGA Data Type Codes:
			See: http://www.paulbourke.net/dataformats/tga/
			0  -  No image data included.
			1  -  Uncompressed, color-mapped images.
			2  -  Uncompressed, RGB images.
			3  -  Uncompressed, black and white images.
			9  -  Runlength encoded color-mapped images.
		   10  -  Runlength encoded RGB images.
		   11  -  Compressed, black and white images.
		   32  -  Compressed color-mapped data, using Huffman, Delta, and
				  runlength encoding.
		   33  -  Compressed color-mapped data, using Huffman, Delta, and
				  runlength encoding.  4-pass quadtree-type process.	*/

		// Get TGA format header.
		unsigned char TGAheader[12];
		fread(TGAheader, sizeof(unsigned char), 12, f);
		char idlength = TGAheader[0];
		char colourmaptype = TGAheader[1];
		char datatypecode = TGAheader[2];
		
		// This function can only handle image type 2 and 10 for now.
		if (datatypecode != 2 && datatypecode != 10)
			return;	

		int colourmaporigin = (TGAheader[4] << 8) ^ TGAheader[3];
		int colourmaplength = (TGAheader[6] << 8) ^ TGAheader[5];
		char colourmapdepth = TGAheader[7];
		int x_origin = (TGAheader[9] << 8) ^ TGAheader[8];
		int y_origin = (TGAheader[11] << 8) ^ TGAheader[10];

		// Get info format header.
		unsigned char header[6];
		fread(header, sizeof(unsigned char), 6, f);
		*cols = (header[1] << 8) ^ header[0];
		*rows = (header[3] << 8) ^ header[2];
		int colorDepth = header[4];
		int bytesPerPixel = (colorDepth / 8);
		int imagedescriptor = header[5];

		// Read the image. This code borrowed from Paul Bourke at:
		// http://www.paulbourke.net/dataformats/tga/tgatest.c.
		int bytes2read = bytesPerPixel;
		unsigned char p[5];
		int n = 0;
		int i = 0;
		int j = 0;
		std::vector<std::array<unsigned char, 4>> temp(*rows * *cols);
		while (n < (*cols) * (*rows)) {
			if (datatypecode == 2) {/* Uncompressed */
				if (fread(p, 1, bytes2read, f) != bytes2read) {
					fprintf(stderr, "Unexpected end of file at pixel %d\n", i);
					return;
				}
				MergeBytes(&temp[n], p, bytes2read);
				n++;
			}
			else if (datatypecode == 10) {/* Compressed */
				if (fread(p, 1, bytes2read + 1, f) != bytes2read + 1) {
					fprintf(stderr, "Unexpected end of file at pixel %d\n", i);
					return;
				}
				j = p[0] & 0x7f;
				MergeBytes(&(temp[n]), &(p[1]), bytes2read);
				n++;
				if (p[0] & 0x80) {/* RLE chunk */
					for (i = 0; i < j; i++) {
						MergeBytes(&(temp[n]), &(p[1]), bytes2read);
						n++;
					}
				}
				else { /* Normal chunk */
					for (i = 0; i < j; i++) {
						if (fread(p, 1, bytes2read, f) != bytes2read) {
							fprintf(stderr, "Unexpected end of file at pixel %d\n", i);
							return;
						}
						MergeBytes(&(temp[n]), p, bytes2read);
						n++;
					}
				}
			}
		}
		fclose(f);

		// Copy pixel data to array.
		for (int idx = 0; idx < temp.size(); idx++)
		{
			if (bytesPerPixel == 1)
				pixels->push_back(static_cast<real>(temp[i][0]));
			else if (bytesPerPixel == 3 || bytesPerPixel == 4)
			{
				std::vector<real> vec = std::vector<real>
				{ 
					static_cast<real>(temp[idx][0]),
					static_cast<real>(temp[idx][1]),
					static_cast<real>(temp[idx][2])
				};
				pixels->push_back(convertRGBtoGray(vec));
			}
		}

		temp.clear();
		fclose(f);
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
					Vector rgb = convertGrayToRGB(val);
					fwrite(rgb.get()->data(), sizeof(unsigned char), 3, f);
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
		std::ofstream ofile(filename, std::ios::out | std::ios::binary);
		if (!ofile.is_open()) 
			return;		

		ofile << "P6" << "\n";
		ofile << cols << " ";
		ofile << rows << "\n";
		ofile << "255" << "\n";
		void* img = (void*)pixels->data(); 
		unsigned char* temp = new unsigned char[pixels->size()*3];
		for (int i = 0; i < pixels->size(); i++) 
		{
			// Convert gray->RGB, then swap bytes from RGB->BGR.
			Vector val = convertGrayToRGB((*pixels)[i]);
			temp[i*3+0] = static_cast<unsigned char>(val[0]);
			temp[i*3+1] = static_cast<unsigned char>(val[1]);
			temp[i*3+2] = static_cast<unsigned char>(val[2]);
		}
		ofile.write((char*)temp, pixels->size()*3);
		delete[] temp;
		ofile.close();
	}

	void saveTGA(std::string filename, int rows, int cols, int colorDepth,
		std::vector<real>* pixels)
	{
		FILE* f = fopen(filename.c_str(), "wb");

		// Transfer elements to pixel array.
		std::vector<unsigned char> temp(pixels->size() * 3, 0);
		for (int i = 0; i < pixels->size(); i++)
		{
			// Convert back from Gray to RGB and write pixels to file.
			Vector val = convertGrayToRGB((*pixels)[i]);
			temp[i * 3 + 0] = static_cast<unsigned char>(val[0]);
			temp[i * 3 + 1] = static_cast<unsigned char>(val[1]);
			temp[i * 3 + 2] = static_cast<unsigned char>(val[2]);
		}

		// Make headers, write all to file.
		unsigned char TGAheader[12] = { 0,0,2,0,0,0,0,0,0,0,0,0 };
		unsigned char header[6] = { cols % 256, cols / 256,
			rows % 256, rows / 256, colorDepth, 0 };
		fwrite(TGAheader, sizeof(unsigned char), 12, f);
		fwrite(header, sizeof(unsigned char), 6, f);
		fwrite(temp.data(), sizeof(unsigned char), temp.size(), f);
		fclose(f);
	}
}