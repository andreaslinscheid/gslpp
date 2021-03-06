#!/bin/bash
#
#	This script generates the 'tests.cpp' file that exectues all tests defined in
#	the <modules>/test/ folder. The format to register a test is given below.
#
#
#	To include a test into the test suite do the following:
#	 - In a new module that has no 'test/Test.h' file:
#		* create the folder 'test/'
#	 	* write the class with the header 'Test.h' that has a public routine:
#			run_test()
#		  This routine will be executed and should call all tests in the module.
#	 - In a module that has the 'test/Test.h' file according to the above:
#		* from the run_test() implementation call the test in whatever way that is
#		  convinient.
#

#
#	find all the 'Test.h' files and prepare their inclusion
#
newline=$'\n'
tab=$'\t'
testFolders=`find ../gslpp/ -type d | grep -P "test$"`
includes=""
executionLines=""
for i  in $testFolders; do
	moduleRelPath=${i%%"/test"}
	module=${moduleRelPath##"../gslpp/"}
	if [[ -f "$i/Test.h" ]]; then
		include=`echo -e "#include \"gslpp/$module/test/Test.h\""`
		declaretion=`echo -e "gslpp::$module::RunTest $module""_test;"`	
		execution=`echo -e "$module""_test.run_test();"`
		#
		#	add the lines to the respective collection
		#
		includes=${includes}${newline}${include}
		executionLines=${executionLines}${newline}${newline}${tab}"//Test of the module $module"
		executionLines=${executionLines}${newline}${tab}${declaretion}${newline}${tab}${execution}
	else
		echo -e "\n WARNING:"
		echo -e "Could not find tests for module $module\n"
	fi
done
#
#	generate the final cpp file
#
timestamp=`date`
cat > tests.cpp << EOF
/*
 * tests.cpp
 *  
 *  Tests all module of the gslpp library that have the 'test/' folder
 *
 *  Generated by generate_tests.sh
 *  DO NOT MODIFY - ALL CHANGES WILL BE LOST
 *  
 *  Time where this file was generated is $timestamp
 */

#include <iostream>
$includes

int main (int argc, char *argv[]){
	std::cout << "Starting test of library : \\n\\n" << std::endl;
	std::cout << "\\tThe test suit was generated at $timestamp" << std::endl;
	$executionLines
	
	
	std::cout << "\\n\\nTest finished" << std::endl;
	return 0;
}
EOF
