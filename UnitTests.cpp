#include "UnitTests.h"
#include "MathParser.h"
#include "StringUtils.h"
#include "Parsing.h"
#include <iostream>
#include <ctime>

namespace MathParser
{
	bool UnitTests::replaceStringTest()
	{
		// Set up test variables.
		srand(clock());
		std::string testStr;
		int charsPerString = 128;
		testStr.resize(charsPerString);
		unsigned int number_of_tests = 100;

		// Set ascii char range for testing.
		char maxChar = '}';
		char minChar = '!';
		char extent = maxChar - minChar;

		for (int n = 0; n < number_of_tests; ++n)
		{
			// Set up variables for replacement.			
			char toReplaceChar = (rand() % extent) + minChar;
			std::string toReplace = " ";
			toReplace[0] = toReplaceChar;

			char replaceWithChar = (rand() % extent) + minChar;
			while (replaceWithChar == toReplaceChar)// Ensure toReplace != replaceWith
				replaceWithChar = (rand() % extent) + minChar;
			std::string replaceWith = " ";
			replaceWith[0] = replaceWithChar;

			// Set up random string.
			for (int i = 0; i < charsPerString; ++i)
			{
				testStr[i] = (rand() % extent) + minChar;
			}

			// Do replacement.
			testStr = replaceString(testStr, toReplace, replaceWith);

			// Check results.
			if (testStr.find(toReplace) != std::string::npos || testStr.length() != charsPerString)
			{
				std::cout << "starting string: " << testStr << std::endl;
				std::cout << "to replace: " << toReplace << std::endl;
				std::cout << "replace with: " << replaceWith << std::endl;
				return false;
			}
		}

		return true;
	}

	bool UnitTests::countOperatorsTest()
	{
		std::string testStr = "(0.2-1)";
		int idx = -1;
		if ((idx = countOperators(testStr)) != 1)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		testStr = "(0.2+-1)";
		if ((idx = countOperators(testStr)) != 1)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		testStr = "2^(0.2-1)";
		if ((idx = countOperators(testStr)) != 2)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::countComplexNumbersTest()
	{
		std::string testStr = "(0.2-1i)";
		int idx = -1;
		if ((idx = countComplexNumbers(testStr)) != 1)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		testStr = "(2+i) * 23-2i";
		if ((idx = countComplexNumbers(testStr)) != 2)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 2" << std::endl;
			return false;
		}
	}

	bool UnitTests::findNearestRealTest()
	{
		std::string testStr = "(24/2) + 2";
		int idx = -1;
		if ( (idx=findNearestReal(testStr)) != 1)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		testStr = " . --4  ";
		if ((idx = findNearestReal(testStr)) != 2)
		{
			std::cout << "test result: " << idx << std::endl;
			std::cout << "expected result: 2" << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::getNearestRealTest()
	{
		std::string testStr = "(24/2) + 2";
		if ( (testStr = getNearestReal(testStr)) != "24")
		{
			std::cout << "test result: " << testStr << std::endl;
			std::cout << "expected result: 1" << std::endl;
			return false;
		}

		testStr = " . --4  ";
		if ( (testStr = getNearestReal(testStr)) != "-4")
		{
			std::cout << "test result: " << testStr << std::endl;
			std::cout << "expected result: 4" << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::parseNInputNumbersTest()
	{
		real acceptableError = 0.001;
		unsigned int input_numbers = 50;
		unsigned int number_of_tests = 500;

		// Test 1: apply arithmetic to two random numbers.
		for (int i = 0; i < number_of_tests; ++i)
		{
			// Create test string of numbers
			std::string testStr;
			std::vector<real> testVals;
			std::vector<real> results;
			for (int j = 0; j < input_numbers; ++j)
			{
				// Create random test values.
				real testVal = real(rand()) / 9999.999;
				testVals.push_back(testVal);
				testStr.append(std::to_string(testVal));
				if(j < input_numbers-1)
					testStr.append(",");
			}

			// Parse numbers.
			results = ParseNInputNumbers(testStr);

			// Check results first for vector size, then for term-wise parsing error.
			bool isCorrect = results.size() == testVals.size();
			if (!isCorrect)
			{
				std::cout << "TERM NUMBER MISMATCH" << std::endl;
				std::cout << "resultant terms: " << results.size() << std::endl;
				std::cout << "expected terms: " << input_numbers << std::endl;

				// Find mismatching terms to locate error.
				int terms = (results.size() < testVals.size()) ? results.size() : testVals.size();
				for (int j = 0; j < terms; ++j)
				{
					if (std::abs(results[j] - testVals[j]) > acceptableError)
					{
						std::cout << "Term mismatch at term #" << j << std::endl;
						std::cout << results[j] << "," << testVals[j] << std::endl;
						break;
					}
				}
			}

			for (int j = 0; isCorrect && (j < input_numbers); ++j)
			{
				if (std::abs(results[j] - testVals[j]) > acceptableError)
				{
					isCorrect = false;
					std::cout << "TERM PRECISION MISMATCH" << std::endl;
					std::cout << "test value: " << results[j] << std::endl;
					std::cout << "expected value: " << testVals[j] << std::endl;
					break;
				}
			}

			// If results are false, end tests early.
			if (!isCorrect)			
				return false;			
		}

		// Test 2: try random junk strings with no numbers.
		std::string testStr = "(,)";
		std::vector<real> results = ParseNInputNumbers(testStr);
		if (results.size() > 0)
		{
			std::cout << "test values: " << results.size() << std::endl;
			std::cout << "expected values: 0" << std::endl;
			return false;
		}

		return true;
	}


	bool UnitTests::calculateArithmeticTest()
	{
		real allowableError = 0.0001;
		unsigned int number_of_tests = 500;
		char operations[4] = { '+', '-', '*', '/' };

		// Test 1: apply arithmetic to two random numbers.
		for (int i = 0; i < number_of_tests; ++i)
		{
			std::string operationString = " ";
			std::string testStr;
			real testVal1, testVal2;
			real testResult = 0;
			real expectedResult = 0;

			// Create two random test values.
			testVal1 = real(rand()) / 9999.999;
			testVal2 = real(rand()) / 9999.999;
			char operation = operations[rand() % 4];
			operationString[0] = operation;

			// Create the test string.
			testStr.append(std::to_string(testVal1));
			while (rand() % 2)
				testStr.append(" ");
			testStr.append(operationString);
			testStr.append(std::to_string(testVal2));
			while (rand() % 2)
				testStr.append(" ");

			// Do arithmetic.
			testResult = calculateArithmetic(testStr);

			// Get actual result.
			switch (operation)
			{
			case '+':
				expectedResult = testVal1 + testVal2;
				break;
			case '-':
				expectedResult = testVal1 - testVal2;
				break;
			case '*':
				expectedResult = testVal1 * testVal2;
				break;
			case '/':
				expectedResult = testVal1 / testVal2;
				break;
			}

			// Check results.
			if (std::abs(testResult - expectedResult) > allowableError)
			{
				std::cout << "test string: " << testStr << std::endl;
				std::cout << "expected result: " << expectedResult << std::endl;
				std::cout << "actual result: " << testResult << std::endl;
				return false;
			}
		}

		return true;
	}

	bool UnitTests::processTwoTermOperationTest()
	{
		std::string testStr = "(-1.25 * 4)";
		std::vector<std::string> expr = loadExpressionsIntoArray(testStr);
		real result = ParseInputNumber(processTwoTermOperation(testStr,expr[0]));
		if (result != -5)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result: -5" << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::parserTest()
	{
		MathParser parser;
		parser.precision = 30;
		real acceptableError = 0.0001;

		// Test 1: simple arithmetic.
		std::string testStr = "(2 + 2) - 4";
		real result = ParseInputNumber(parser.parse(testStr));
		if (result != 0)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result: 0" << std::endl;
			return false;
		}

		// Test 2: order of operations test.
		testStr = "4 + (-1.25 * 4) / 2 + 0.2";
		result = ParseInputNumber(parser.parse(testStr));
		if (std::abs(result - 1.7) > acceptableError)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result: 1.3" << std::endl;
			return false;
		}

		// Test 3:
		testStr = "(24 / 2) * 0.5";
		result = ParseInputNumber(parser.parse(testStr));
		if (std::abs(result - 6) > acceptableError)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result: 6" << std::endl;
			return false;
		}

		// Test 4:
		testStr = "(24/2)*0.5";
		result = ParseInputNumber(parser.parse(testStr));
		if (std::abs(result - 6) > acceptableError)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result: 6" << std::endl;
			return false;
		}

		// Test 5:
		testStr = "(13 * 2) / 3";
		result = ParseInputNumber(parser.parse(testStr));
		if (std::abs(result - 8.6666666667) > acceptableError)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result:  8.6666666667" << std::endl;
			return false;
		}

		// Test 6:
		testStr = "2^(0.2-1)";
		result = ParseInputNumber(parser.parse(testStr));
		if (std::abs(result - 0.5743491775) > acceptableError)
		{
			std::cout << "test result: " << result << std::endl;
			std::cout << "expected result:  0.5743491775" << std::endl;
			return false;
		}

		//// Test 7: Parse arithmetic with complex numbers
		//testStr = "(2+3i)*(1+1i)";
		//complex resultComplex = ParseComplexNumber(parser.parse(testStr));
		//if (std::abs((resultComplex - std::complex<real>(-1.0,6.0)).real()) > acceptableError
		// && std::abs((resultComplex - std::complex<real>(-1.0,6.0)).imag()) > acceptableError)
		//{
		//	std::cout << "test result: " << resultComplex << std::endl;
		//	std::cout << "expected result: -1+6i" << std::endl;
		//	return false;
		//}

		return true;
	}

	void UnitTests::runAll()
	{
		bool passed = false;

		// Run replaceString test.
		if (!(passed = replaceStringTest()))
		{
			std::cout << "ERROR! Replace string test failed." << std::endl;
			return;
		}

		// Find nearest real test.
		if (!(passed = findNearestRealTest()))
		{
			std::cout << "ERROR! Find nearest real test failed." << std::endl;
			return;
		}

		// Get nearest real test.
		if (!(passed = getNearestRealTest()))
		{
			std::cout << "ERROR! Get nearest real test failed." << std::endl;
			return;
		}

		// Find count operators test.
		if (!(passed = countOperatorsTest()))
		{
			std::cout << "ERROR! count operators test failed." << std::endl;
			return;
		}

		// Find count complex numbers test.
		if (!(passed = countComplexNumbersTest()))
		{
			std::cout << "ERROR! count complex numbers test failed." << std::endl;
			return;
		}

		// Run parseNInputNumbers test.
		if (!(passed = parseNInputNumbersTest()))
		{
			std::cout << "ERROR! parseNInputNumbers test failed." << std::endl;
			return;
		}

		// Run processTwoTermOperation test.
		if (!(passed = processTwoTermOperationTest()))
		{
			std::cout << "ERROR! processTwoTermOperation test failed." << std::endl;
			return;
		}

		// Run calculateArithmetic test.
		if (!(passed = calculateArithmeticTest()))
		{
			std::cout << "ERROR! Calculate arithmetic test failed." << std::endl;
			return;
		}

		// Run parser test.
		if (!(passed = parserTest()))
		{
			std::cout << "ERROR! Parser test failed." << std::endl;
			return;
		}

		std::cout << "All tests successfully passed." << std::endl;
	}
}