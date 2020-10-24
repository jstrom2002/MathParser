#include "MathParser.h"
#include "UnitTests.h"
#include <iostream>

int main(int argc, char** argv)
{
	// If argv[1] = 'TEST', run program in test mode.
	if (argc > 1 && std::string(argv[1]) == "TEST")
	{
		MathParser::UnitTests tests;
		tests.runAll();
	}
	else
	{
		bool quitProgram = false;
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "//           MATH PARSER          //" << std::endl;
		std::cout << "//                                //" << std::endl;
		std::cout << "// How to use: type any arith-    //" << std::endl;
		std::cout << "// metic calculation in or enter  //" << std::endl;
		std::cout << "// a function as 'f(x) = [terms]' //" << std::endl;
		std::cout << "// Enter polynomials similarly    //" << std::endl;
		std::cout << "// as 'p(x) = [terms]'            //" << std::endl;
		std::cout << "//--------------------------------//" << std::endl;
		std::cout << "    <enter 'q' to quit program>     " << std::endl;
		std::cout << std::endl;
		char inputTxt[256];
		MathParser::MathParser parser;
		std::string stringOutput = "";

		while (!quitProgram)
		{
			// Print functions saved in memory.
			std::string functionStr = parser.printSavedFunctions();
			if (functionStr.length())
			{
				std::cout << std::endl;
				std::cout << "Saved Functions:" << std::endl;
				std::cout << functionStr;
				std::cout << "================" << std::endl;
			}

			// Print polynomials saved in memory.
			std::string polyStr = parser.printSavedPolynomials();
			if (polyStr.length())
			{
				std::cout << std::endl;
				std::cout << "Saved Polynomials:" << std::endl;
				std::cout << polyStr;
				std::cout << "================" << std::endl;
			}

			std::cout << "Enter equation: ";
			std::cin.getline(inputTxt, 256);
			if (inputTxt[0] == 'q' && inputTxt[1] == '\0')
			{
				quitProgram = true;
				return 0;
			}
			else if( (stringOutput = parser.parse(inputTxt)).length() )
			{
				std::cout << "result: " << stringOutput << std::endl;
			}
		}
	}
}