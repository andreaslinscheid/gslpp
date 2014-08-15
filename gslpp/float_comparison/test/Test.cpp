/*
 * Test.cpp
 *
 *  Created on: Jun 8, 2014
 *      Author: alinsch
 *
 *	A small program to test the error handling object
 */
#include <iostream>
#include "gslpp/float_comparison/test/Test.h"
#include "gslpp/float_comparison/FloatComparison.h"

namespace gslpp{
namespace float_comparison{

void RunTest::run_test(){
	std::cout << "\n\nTest of the Float comparison methods" <<std::endl;

	//equal numbers compare equal
	double a = 1e-5;
	double b = 1e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,4) )
		std::cout << "Test failure: equal number do not agree in first 4 digits" <<std::endl;

	//equal number still are equal
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,0) )
		std::cout << "Test failure: equal number do not agree in first 0 digits" <<std::endl;

	//number that differ in the first digit
	a = 1e-5;
	b = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,0) )
		std::cout << "Test failure: a and be differ in first digit so agreement to 0 digits is fine" <<std::endl;
	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,1) )
		std::cout << "Test failure: a and be differ in first digit" <<std::endl;

	//test the threshold functionality
	a = 1e-5;
	b = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(a,b,4,1e-4) )
		std::cout << "Test failure: threshold does not remove the comparison" <<std::endl;

	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(a,b,2,1.5e-5) )
		std::cout << "Test failure: threshold not does not remove the comparison if one number is larger" <<std::endl;

	a = 0.0 / 0.0;//NaN
	if ( not gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: NaN check does not work" <<std::endl;

	a = - 0.0 / 0.0;//NaN
	if ( not gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: NaN check does not work" <<std::endl;

	a = 1.0 / 0.0;//infinity
	if ( gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: Infinity is not NaN" <<std::endl;

	std::cout << "Test success" <<std::endl;
}

};
};
