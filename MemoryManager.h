#pragma once
#include "Types.h"
#include "ComplexNumber.h"
#include "Function.h"
#include "Polynomial.h"
#include "Matrix.h"
#include "Vector.h"


namespace MathParser
{
	class MemoryManager;
	class UserInterface;

	class MemoryManager
	{
	public:
		MemoryManager();

		void add(Function object);
		void add(Polynomial object);
		void add(Matrix object);
		void add(Vector object);
		void clear();
		void loadFile(std::string str);
		void print(int precision);
		void saveFile(std::string str);
		size_t size();

		friend class MathParser;
		friend class UserInterface;

	private:
		real LAST;
		std::vector<Function> functions;
		std::vector<Matrix> matrices;
		std::vector<Polynomial> polynomials;
		std::vector<Vector> vectors;

		std::string printSavedFunctions(int precision);
		std::string printSavedMatrices(int precision);
		std::string printSavedPolynomials(int precision);
		std::string printSavedVectors(int precision);
	};
}