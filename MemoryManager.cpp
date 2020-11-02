#include "MemoryManager.h"
#include "StringUtils.h"
#include "FileIO.h"
#include <iostream>

namespace MathParser
{
	MemoryManager::MemoryManager() : LAST(0)
	{
	}

	void MemoryManager::add(Function object){functions.push_back(object);}
	void MemoryManager::add(Matrix object){matrices.push_back(object);}
	void MemoryManager::add(Polynomial object){polynomials.push_back(object);}
	void MemoryManager::add(Vector object){vectors.push_back(object);}

	void MemoryManager::clear()
	{
		functions.clear();
		matrices.clear();
		polynomials.clear();
		vectors.clear();
	}

	void MemoryManager::print(int precision)
	{
		if (!size())
			return;

		std::cout << "Saved:" << std::endl;
		std::cout << printSavedFunctions(precision);
		std::cout << printSavedMatrices(precision);
		std::cout << printSavedPolynomials(precision);
		std::cout << printSavedVectors(precision);
		std::cout << std::endl;
	}

	std::string MemoryManager::printSavedFunctions(int precision)
	{
		std::string str = "";
		if (!functions.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < functions.size(); ++i)
		{
			std::string funcName = "f" + std::to_string(i) + " - Function";
			str.append(funcName + newline);
		}

		return str;
	}

	std::string MemoryManager::printSavedMatrices(int precision)
	{
		std::string str = "";
		if (!matrices.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < matrices.size(); ++i)
		{
			std::string funcName = "M" + std::to_string(i) + " - Matrix size [" +
				std::to_string(matrices[i].rows) + "," + 
				std::to_string(matrices[i].columns) + "]";
			str.append(funcName + newline);
		}

		return str;
	}

	std::string MemoryManager::printSavedPolynomials(int precision)
	{
		std::string str = "";
		if (!polynomials.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < polynomials.size(); ++i)
		{
			std::string funcName = "p" + std::to_string(i) + " - Polynomial size "
				+ std::to_string(polynomials[i].size());
			str.append(funcName + newline);
		}

		return str;
	}

	std::string MemoryManager::printSavedVectors(int precision)
	{
		std::string str = "";
		if (!vectors.size())
			return str;

		std::string newline = "\n";
		for (int i = 0; i < vectors.size(); ++i)
		{
			std::string funcName = "V" + std::to_string(i) + " - Vector size " + 
				std::to_string(vectors.size());
			str.append(funcName + newline);
		}

		return str;
	}

	void MemoryManager::saveFile(std::string str)
	{
	}

	void MemoryManager::loadFile(std::string str)
	{
		// Remove ')' around outside of filename.
		while (str[str.length() - 1] == ')')
			str = str.substr(0, str.length() - 1);

		// Remove all chars before "(" from from beginning of filename.
		str = str.substr(str.find("(") + 1);

		// Handle file loading according to file extension.
		std::string ext = getExtension(str);
		if (str.find(".csv") != std::string::npos)
			matrices.push_back(csvread(str));
		else if (str.find(".bmp") != std::string::npos || 
			str.find(".ppm") != std::string::npos || 
			str.find(".tga") != std::string::npos)
				matrices.push_back(imread(str));
	}

	size_t MemoryManager::size()
	{
		return functions.size() + polynomials.size() + matrices.size() + 
			vectors.size();
	}
}