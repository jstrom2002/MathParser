#include "UI.h"
#include "Parsing.h"
#include "UnitTests.h"
#include <string>
#include <iostream>

namespace MathParser
{
	void UserInterface::printSavedMemory()
	{
		// Print functions saved in memory.
		functionString = parser.printSavedFunctions();
		if (functionString.length())
		{
			std::cout << std::endl;
			std::cout << "Saved Functions:" << std::endl;
			std::cout << functionString;
			std::cout << "================" << std::endl;
		}

		// Print polynomials saved in memory.
		polynomialString = parser.printSavedPolynomials();
		if (polynomialString.length())
		{
			std::cout << std::endl;
			std::cout << "Saved Polynomials:" << std::endl;
			std::cout << polynomialString;
			std::cout << "================" << std::endl;
		}
	}

	void UserInterface::interpretInput()
	{
		// System commands:
		if (!strcmp(inputText, "exit"))
			windowClose = true;
		else if (!strcmp(inputText, "clear"))
			parser.clear();
		else if (!strcmp(inputText, "help"))
			printHelp();
		else if (!strcmp(inputText, "test"))
		{
			UnitTests tests;
			tests.runAll();
		}
		else if (std::string(inputText).find("precision") != std::string::npos)
		{
			std::string str(inputText);
			str = str.substr(str.find("precision") + 1);
			str = str.substr(str.find("=") + 1);
			parser.precision = ParseInputNumber(str);
		}

		// Else parse mathematical statements.
		else if ((outputText = parser.parse(inputText)).length())		
			std::cout << "result: " << outputText << std::endl;		
	}

	void UserInterface::printHelp()
	{
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "//           MATH PARSER          //" << std::endl;
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "List of system commands:" << std::endl;
		std::cout << "clear" << std::endl;
		std::cout << "exit" << std::endl;
		std::cout << "help" << std::endl;
		std::cout << "test" << std::endl;
		std::cout << std::endl;

		std::cout << "List of user-defined functions:" << std::endl;
		for (int i = 0; i < functionList.size(); ++i)
			std::cout << functionList[i] << std::endl;
		std::cout << std::endl;

		std::cout << "To use last entered value in an expression, type 'LAST'." << std::endl;
		std::cout << "For example, a user could type 'LAST + 2'" << std::endl;
		std::cout << "To change precision, type 'precision = n'" << std::endl;
		std::cout << "To declare a function, type 'f(x,y,..) = [equation terms]'" << std::endl;
		std::cout << "To declare a polynomial, type 'p(x) = [polynomial terms like 'x^2 + x']'" << std::endl;
		std::cout << "To reference a saved mathematical object (like a function), use the first letter of the"
			"object's name, followed by the index number. For example, to evaluate function 1 at x = 2.3, type"
			"f1.evaluate(2.3)" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for function class: " << std::endl;
		std::cout << ".evaluate(real)" << std::endl;
		std::cout << ".derivative(real D, real x, real epsilon)" << std::endl;
		std::cout << ".integral(real a, real b, int abscissaPoints)" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for polynomial class: " << std::endl;
		std::cout << ".evaluate(real)" << std::endl;
		std::cout << ".realRoots()" << std::endl;
		std::cout << ".roots()" << std::endl;
		std::cout << std::endl;
	}

	void UserInterface::run()
	{
		std::cout << "<enter 'help' for list of commands>" << std::endl;
		std::cout << std::endl;	

		while (!windowClose)
		{
			printSavedMemory();
			std::cout << "Enter equation: ";
			std::cin.getline(inputText, 256);			
			interpretInput();
		}
	}
}