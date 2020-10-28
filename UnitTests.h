/**
*  Class for running unit tests for parsing functions.
*/

#pragma once

namespace MathParser
{
	class UnitTests 
	{
	public:
		void runAll();

	private:
		bool calculateArithmeticTest();
		bool countOperatorsTest();
		bool countComplexNumbersTest();
		bool findNearestRealTest();
		bool getNearestRealTest();
		bool imageProcessingTest();
		bool matrixClassTest();
		bool operatorTest();
		bool parseNInputNumbersTest();
		bool parserTest();
		bool polynomialArithmeticTest();
		bool processTwoTermOperationTest();
		bool replaceStringTest();
		bool rootTest();
	};
}