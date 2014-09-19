/*
 * Test.cpp
 *
 *  Created on: Jun 8, 2014
 *      Author: alinsch
 *
 *	A small program to test the error handling object
 */
#include <iostream>
#include "gslpp/data_interpolation/CubeSpline.h"
#include "gslpp/data_interpolation/CubicPolynomial.h"
#include "gslpp/data_interpolation/test/Test.h"

namespace gslpp{
namespace data_interpolation{

void RunTest::run_test(){
	_allSuccess = true;
	std::cout << "\n\nStarting tests of the data_interpolation namespace" <<std::endl;

	test_CubicPolynomial<double>();
	test_CubicPolynomial<float>();

	test_CubeSpline<double>();
	test_CubeSpline<float>();

	test_HermitePolynomial<double>();
	test_HermitePolynomial<float>();

	test_MonotoneCubeHermiteSpline<double>();
	test_MonotoneCubeHermiteSpline<float>();

	std::cout << "Test of the data_interpolation namespace " <<
			(_allSuccess ? "successfull!" : "not sucessfull!" ) << "\n\n" <<std::endl;
}

template<>
std::string RunTest::nameOfTypeTrait<double> () const {
	return "double";
}
template<>
std::string RunTest::nameOfTypeTrait<float> () const {
	return "float";
}
template<>
double RunTest::accuracyGoal<double> () const {
	return 1e-12;
}
template<>
float RunTest::accuracyGoal<float> () const {
	return 1e-5;
}

};/* namespace data_interpolation */
};/* namespace gslpp */
