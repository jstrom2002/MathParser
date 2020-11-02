#pragma once
#include "MathParser.h"

namespace MathParser
{
	class MemoryManager;

	class UserInterface
	{
	public:
		UserInterface() : windowClose(false) {}
		void run();

		friend class MemoryManager;

	private:
		bool windowClose;
		char inputText[256];
		std::string outputText;
		std::string functionString;
		std::string matricesString;
		std::string polynomialString;
		MathParser parser;

		void printHelp();
		void printSavedMemory();
		void interpretInput();
	};
}