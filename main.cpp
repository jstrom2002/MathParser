#include "UnitTests.h"
#include "UI.h"

int main(int argc, char** argv)
{
	// If argv[1] = 'TEST', run program in test mode.
	if (argc > 1 && !strcmp(argv[1],"TEST"))
	{
		MathParser::UnitTests tests;
		tests.runAll();
	}
	else
	{
		MathParser::UserInterface ui;
		ui.run();		
	}
	return 0;
}