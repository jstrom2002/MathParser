#include "UI.h"
#include "Parsing.h"
#include "UnitTests.h"
#include <string>
#include <iostream>

namespace MathParser
{
	void UserInterface::interpretInput()
	{		
		std::string str(inputText);

		// System commands:
		if (!strcmp(inputText, "exit") || !strcmp(inputText, "quit"))
			windowClose = true;
		else if (!strcmp(inputText, "clc"))
			parser.clc();
		else if (!strcmp(inputText, "clr"))
		{
			parser.clc();
			parser.memory.clear();
		}
		else if (!strcmp(inputText, "clear"))
			parser.clear();
		else if (str.find("disp(") != std::string::npos)
			outputText = parser.disp(inputText);
		else if (!strcmp(inputText, "help"))
			printHelp();
		else if (str.find("load(") != std::string::npos)
		{
			parser.memory.loadFile(inputText);
		}	
		else if (str.find("save(") != std::string::npos)
		{
			parser.memory.saveFile(inputText);
		}
		else if (str.find("type ") != std::string::npos)
		{
			outputText = parser.type(inputText);
		}
		else if (!strcmp(inputText, "test"))
		{
			UnitTests tests;
			tests.runAll();
		}
		else if (str.find("precision") != std::string::npos)
		{
			std::string str(inputText);
			str = str.substr(str.find("precision") + 1);
			str = str.substr(str.find("=") + 1);
			parser.precision = ParseInputNumber(str);
		}

		// Else parse mathematical statements.
		else if ((outputText = parser.parse(inputText)).length())
		{
		}
			
		std::cout << "result: " << outputText << std::endl;		
		std::cout << std::endl;		
	}

	void UserInterface::printHelp()
	{
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "//           MATH PARSER          //" << std::endl;
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "List of system commands:" << std::endl;
		std::cout << "clc" << std::endl;
		std::cout << "clear" << std::endl;
		std::cout << "clr" << std::endl;
		std::cout << "disp" << std::endl;
		std::cout << "exit" << std::endl;
		std::cout << "load(filename)" << std::endl;
		std::cout << "help" << std::endl;
		std::cout << "quit" << std::endl;
		std::cout << "save(filename)" << std::endl;
		std::cout << "test" << std::endl;
		std::cout << "type" << std::endl;
		std::cout << std::endl;

		std::cout << "List of user-defined functions:" << std::endl;
		for (int i = 0; i < functionList.size(); ++i)
			std::cout << functionList[i] << std::endl;
		std::cout << std::endl;

		std::cout << "To use last entered value in an expression, type 'LAST'." << std::endl;
		std::cout << "For example 'LAST + 2'." << std::endl;
		std::cout << "To change precision, type 'precision = n'" << std::endl;
		std::cout << "To declare a function, type 'f(x,y,..) = [equation terms]'" << std::endl;
		std::cout << "To declare a matrix, type rows ending in semicolons, surrounded by brackets ie '[1,0,0;0,1,0;0,0,1]'" << std::endl;
		std::cout << "To declare a polynomial, type 'p(x) = [polynomial terms like 'x^2 + x']'" << std::endl;
		std::cout << "To reference a saved mathematical object (like a function), use the first letter of the"
			" object's name, followed by the index number. For example, to evaluate function 1 at x = 2.3, type"
			" f1.evaluate(2.3)" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for Function class: " << std::endl;
		std::cout << ".evaluate(real)" << std::endl;
		std::cout << ".derivative(real D, real x, real epsilon)" << std::endl;
		std::cout << ".integral(real a, real b, int abscissaPoints)" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for Polynomial class: " << std::endl;
		std::cout << ".evaluate(real)" << std::endl;
		std::cout << ".realRoots()" << std::endl;
		std::cout << ".roots()" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for Matrix class: " << std::endl;
		std::cout << ".det()" << std::endl;
		std::cout << ".inverse()" << std::endl;
		std::cout << ".reshape(int, int)" << std::endl;
		std::cout << ".submatrix(int i1, int j1, int i2, j2)" << std::endl;
		std::cout << std::endl;

		std::cout << "List of methods for Vector class: " << std::endl;
		std::cout << std::endl;
	}

	void UserInterface::run()
	{
		std::cout << "<enter 'help' for list of commands>" << std::endl;
		std::cout << std::endl;	

		while (!windowClose)
		{
			parser.memory.print(parser.precision);
			std::cout << "Enter a command: ";
			std::cin.getline(inputText, 256);			
			interpretInput();
		}
	}
}