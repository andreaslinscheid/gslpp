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
		std::cout << "Test failure: equal double number do not agree in first 4 digits" <<std::endl;
	float af = 1e-5;
	float bf = 1e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(af,bf,4) )
		std::cout << "Test failure: equal float number do not agree in first 4 digits" <<std::endl;

	//equal number still are equal
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,0) )
		std::cout << "Test failure: equal number do not agree in first 0 digits" <<std::endl;

	//number that differ in the first digit
	a = 1e-5;
	b = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,0) )
		std::cout << "Test failure: a and b differ in first digit so agreement to 0 digits is fine" <<std::endl;
	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal(a,b,1) )
		std::cout << "Test failure: a and b differ in first digit" <<std::endl;
	af = 1e-5;
	bf = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal(af,bf,0) )
		std::cout << "Test failure: af and bf differ in first digit so agreement to 0 digits is fine" <<std::endl;
	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal(af,bf,1) )
		std::cout << "Test failure: af and bf differ in first digit" <<std::endl;

	//test the threshold functionality
	a = 1e-5;
	b = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(a,b,4,1e-4) )
		std::cout << "Test failure: threshold does not remove the comparison of doubles" <<std::endl;
	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(a,b,2,1.5e-5) )
		std::cout << "Test failure: threshold not does not remove the comparison if one double number is larger" <<std::endl;
	af = 1e-5;
	bf = 2e-5;
	if ( not gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(af,bf,4,1e-4f) )
		std::cout << "Test failure: threshold does not remove the comparison of floats" <<std::endl;
	if ( gslpp::float_comparison::equal_up_to_significant_digits_decimal_above_threshold(af,bf,2,1.5e-5f) )
		std::cout << "Test failure: threshold not does not remove the comparison if one float number is larger" <<std::endl;

	//Test the NaN check functionality
	a = 0.0 / 0.0;//double NaN
	if ( not gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: NaN check does not work for double" <<std::endl;

	af =  0.0 / 0.0;//float NaN
	if ( not gslpp::float_comparison::is_NaN(af) )
		std::cout << "Test failure: NaN check does not work for float" <<std::endl;

	std::complex<float> ac( 1.0 , 0.0 / 0.0);
	if ( not gslpp::float_comparison::is_NaN(ac) )
		std::cout << "Test failure: NaN check does not work for complex float with NaN imag part" <<std::endl;
	ac = std::complex<float>( 0.0 / 0.0 ,  1.0 );
	if ( not gslpp::float_comparison::is_NaN(ac) )
		std::cout << "Test failure: NaN check does not work for complex float with NaN real part" <<std::endl;

	a = - 0.0 / 0.0;//NaN
	if ( not gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: NaN check does not work for double" <<std::endl;

	a = 1.0 / 0.0;//infinity
	if ( gslpp::float_comparison::is_NaN(a) )
		std::cout << "Test failure: Infinity is not NaN" <<std::endl;


	std::cout << "Test success" <<std::endl;
}

};
};
