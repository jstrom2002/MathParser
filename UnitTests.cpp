#include "UnitTests.h"
#include "MathParser.h"
#include "StringUtils.h"
#include "Parsing.h"
#include "Polynomial.h"
#include "Matrix.h"
#include "ImageProcessing.h"
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

	bool UnitTests::polynomialArithmeticTest()
	{
		// Initialize test polynomials.
		Polynomial p1(std::vector<real>{1, 2, 2});//x^2 + 2x + 2.
		Polynomial p2(std::vector<real>{1, 0, -1});//x^2 - 1

		Polynomial p3 = p1 + 2.0;
		std::string testStr = p3.to_string();
		if (testStr != "x^2 + 2x + 4")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p1 + 2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = 2.0 + p1;
		testStr = p3.to_string();
		if (testStr != "x^2 + 2x + 4")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p1 * 2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = p1 + p2;
		testStr = p3.to_string();
		if (testStr != "2x^2 + 2x + 1")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p2 = " << p2.to_string() << std::endl;
			std::cout << "p1 + p2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = p1 - 2.0;
		testStr = p3.to_string();
		if (testStr != "x^2 + 2x")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p1 - 2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = 2.0 - p1;
		testStr = p3.to_string();
		if (testStr != "-x^2 - 2x")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "2 - p1 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = p1 - p2;
		testStr = p3.to_string();
		if (testStr != "2x + 3")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p2 = " << p2.to_string() << std::endl;
			std::cout << "p1 - p2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = p1 * 2.0;
		testStr = p3.to_string();
		if (testStr != "2x^2 + 4x + 4")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p1 * 2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = 3.0 * p1;
		testStr = p3.to_string();
		if (testStr != "3x^2 + 6x + 6")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "3.0 * p1 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = (p1 * p2);
		testStr = p3.to_string();
		if (testStr != "x^4 + 2x^3 + x^2 - 2x - 2")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p2 = " << p2.to_string() << std::endl;
			std::cout << "p1 * p2 = " << p3.to_string() << std::endl;
			return false;
		}

		p3 = (p1 * p2);
		testStr = p3.to_string();
		if (testStr != "x^4 + 2x^3 + x^2 - 2x - 2")
		{
			std::cout << "p1 = " << p1.to_string() << std::endl;
			std::cout << "p2 = " << p2.to_string() << std::endl;
			std::cout << "p1 * p2 = " << p3.to_string() << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::rootTest()
	{	
		// Don't forget the input real array must be from lowest-to-highest
		// exponent, the following initializes the polynomial 'x^2 - 5x + 6'.
		std::vector<real> roots = Polynomial(std::vector<real>{1, -5, 6}).realRoots();
		if (roots.size() != 2 || roots[0] != 3 || roots[1] != 2)
		{
			std::cout << "expected roots: {3,2}" << std::endl;
			std::cout << "test results = {";
			for (int i = 0; i < roots.size(); ++i)
			{
				std::cout << roots[i];
				if(i<roots.size()-1)
					std::cout << ",";
				return false;
			}
			std::cout << "}" << std::endl;
		}

		//roots = Polynomial(std::vector<real>{-2,-1,-2,1}).realRoots();
		//if (roots.size() != 3 || roots[0] != 2 || roots[1] != 1 || roots[2] != -1)
		//{
		//	std::cout << "expected roots: {2,1,-1}" << std::endl;
		//	if (!roots.size())
		//		return false;
		//	std::cout << "test results = {";
		//	for (int i = 0; i < roots.size(); ++i)
		//	{
		//		std::cout << roots[i];
		//		if (i < roots.size() - 1)
		//			std::cout << ",";
		//		return false;
		//	}
		//	std::cout << "}" << std::endl;
		//}

		return true;
	}

	bool UnitTests::operatorTest()
	{
		std::vector<real> p1 = std::vector<real>{
			1,0,-2,
			0,3,-1
		};
		std::vector<real> p2 = std::vector<real>{
			 0, 3,
			-2,-1,
			 0, 4
		};
		std::vector<real> p3 = std::vector<real>{
			0,-5,
		   -6,-7
		};
		std::vector<real> p4 = std::vector<real>{
			3,2,2,
			2,3,-2
		};
		Matrix A1(2,3,p4);
		Matrix A2(2,3,p4);
		
		// Test boolean operators.
		if (A1 != A2)
		{
			std::cout << "Matrix '==' operator failed" << std::endl;
			return false;
		}

		// Test '()' operator. First, check accessing values.
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 3; ++j)
				if (A1(i, j) != p4[i * 3 + j])
				{
					std::cout << "test result: " << A1(i, j) << std::endl;
					std::cout << "expected result: " << p4[i * 3 + j] << std::endl;
					return false;
				}

		// Now check setting values with '()' operator.
		for (int i = 0; i < 2; ++i)
			for (int j = 0; j < 3; ++j)
			{
				A1(i, j) = 0;
				if (A1(i, j))
				{
					std::cout << "test result: " << A1(i, j) << std::endl;
					std::cout << "expected result: 0" << std::endl;
					return false;
				}
			}

		// Check arithmetic operators.
		A1 += A2;
		if (A1 != A2)
		{
			std::cout << "test result: " << A1.to_string() << std::endl;
			std::cout << "expected result: " << A2.to_string() << std::endl;
			return false;
		}

		A1 -= A2;
		Matrix emptyMat = Matrix(A2.rows, A2.columns);
		if (A1 != emptyMat)
		{
			std::cout << "test result: " << A1.to_string() << std::endl;
			std::cout << "expected result: " << emptyMat.to_string() << std::endl;
			return false;
		}

		A2 *= emptyMat;
		if (A2.size())// Check to see if mismatched matrix sizes are caught.
		{
			std::cout << "test result: " << A2.to_string() << std::endl;
			std::cout << "expected result: " << emptyMat.to_string() << std::endl;
			return false;
		}

		// Check results of matrix multiplication of A2 by zero matrix.
		A2 = A1;
		emptyMat = emptyMat.transpose();
		A2 *= emptyMat;
		Matrix A3(A2.rows, A2.columns);
		if (A2 != A3 || A2.max() != 0)
		{
			std::cout << "test result: " << A1.to_string() << std::endl;
			std::cout << "expected result: " << emptyMat.to_string() << std::endl;
			return false;
		}

		// Final matrix multiplication test.
		Matrix M1(2, 3, p1);
		Matrix M2(3, 2, p2);
		Matrix M3(2, 2, p3);
		Matrix Mresult = M1 * M2;
		if (Mresult != M3)
		{
			std::cout << "test result: " << Mresult.to_string() << std::endl;
			std::cout << "expected result: " << M3.to_string() << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::fileIOTest()	
	{
		std::vector<real> p1 = std::vector<real>{
			1, 2, 3,
			4, 5, 6,
			7, 8, 9,
		};
		Matrix M1 = csvread("test.csv");
		Matrix M2(3, 3, p1);

		if (M1 != M2)
		{
			std::cout << "test result: " << M1.to_string() << std::endl;
			std::cout << "expected result: " << M2.to_string() << std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::matrixClassTest()
	{
		std::vector<real> p = std::vector<real>{
			1, 2, 2, 0, 
			0, 1, 3, 1, 
			0, 1, 2, 1, 
			1, 2, 2, -1
		};
		std::vector<real> p2 = std::vector<real>{
			1.0, 0.8, 1.0, 1.0, 0.8, 1.0,
			1.0, 0.5, 0.3, 0.0, 0.5, 1.0,
			1.0, 0.3, 0.2, 0.0, 0.3, 1.0,
			1.0, 0.2, 0.0, 0.0, 0.2, 1.0,
			1.0, 0.3, 0.2, 0.0, 0.3, 1.0,
			1.0, 0.8, 1.0, 1.0, 0.8, 1.0
		};

		std::vector<real> p3 = std::vector<real>{
			1, -2, 2, 4,
			2, -4, 5, 9,
			3, -6, 8, 14,
			5, -10, 12, 22
		};
		std::vector<real> p4 = std::vector<real>{
			3,2,2,
			2,3,-2
		};
		Matrix A(4, 4, p);
		Matrix B(6, 6, p2);
		Matrix C(4, 4, p3);
		Matrix D(2, 3, p4);
		Matrix Mat2 = GaussianKernel2D(5,1);

		// Test: Gauss-Jordan Elimination (can be extremely slow or not converge at all)
		Matrix GJE = A.GaussianElimination();
		Matrix iden = eye(GJE.rows, GJE.columns);
		if (GJE != iden)
		{
			std::cout << "test result: " << GJE.to_string() << std::endl;
			std::cout << "expected result: " << iden.to_string()
				<< std::endl;
			return false;
		}

		return true;
	}

	bool UnitTests::imageProcessingTest()
	{
		// Test 1: load .bmp file, then apply Gaussian blur filter and save results to file.
		Matrix bwImg = imread("cat.bmp");
		Matrix kernel = GaussianKernel2D(5, 1);
		Matrix Mat2 = convolve(kernel, bwImg);
		Mat2.transpose();
		imwrite("test1.bmp", Mat2, 24);
		//BMPtoText("test1.bmp");//also do hex dump of file.

		// Test 2:
		//Matrix Mat = testMatrix();
		//IMAGE bmp2(Mat);

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

		// Run file I/O test.
		if (!(passed = fileIOTest()))
		{
			std::cout << "ERROR! File I/O test failed." << std::endl;
			return;
		}

		// Run operator test.
		if (!(passed = operatorTest()))
		{
			std::cout << "ERROR! Operator test failed." << std::endl;
			return;
		}

		// Run polynomial arithmetic test.
		if (!(passed = polynomialArithmeticTest()))
		{
			std::cout << "ERROR! Polynomial arithmetic test failed." << std::endl;
			return;
		}

		// Run polynomial root test.
		if (!(passed = rootTest()))
		{
			std::cout << "ERROR! Root test failed." << std::endl;
			return;
		}

		// Run matrix class test.
		if (!(passed = matrixClassTest()))
		{
			std::cout << "ERROR! Matrix class test failed." << std::endl;
			return;
		}

		// Run image processing test.
		if (!(passed = imageProcessingTest()))
		{
			std::cout << "ERROR! Image processing test failed." << std::endl;
			return;
		}

		std::cout << "All tests successfully passed." << std::endl;
	}
}