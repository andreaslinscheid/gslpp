/*
 * Test.hpp
 *
 *  Created on: Sep 5, 2014
 *      Author: alinsch
 */
#include "gslpp/data_interpolation/test/Test.h"
#include "gslpp/data_interpolation/CubicPolynomial.h"
#include "gslpp/data_interpolation/CubeSpline.h"

namespace gslpp{
namespace data_interpolation{

template<typename T>
void RunTest::test_CubeSpline() {

	//initialized with the polynomial test data the interpolation must lie on top
	//define some test data
	std::vector<T> xValues;
	std::vector<T> dataSet;
	std::vector<T> derivativeData;
	this->get_polynomial_test_data(xValues,dataSet,derivativeData);

	gslpp::data_interpolation::CubeSpline<T> cubeSpline;
	cubeSpline.initialize(xValues,dataSet);

	gslpp::data_interpolation::CubicPolynomial<T> cpd_testData(xValues.front(),xValues.back(),
			dataSet.front(),dataSet.back(),
			derivativeData.front(),derivativeData.back());

	const size_t numEvals = 100;
	for ( size_t i = 0 ; i < numEvals ; i++){
		T x = cpd_testData.min() + ((cpd_testData.max() - cpd_testData.min() )/numEvals)*i;
		T diffValue = cpd_testData(x) - cubeSpline(x);
		if ( std::fabs(diffValue) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " cubic spline failed: "
					"interpolating test data at step x = " << x << "\n";
			std::cout << "\tData value from spline: " <<  cpd_testData(x) << "; Data value from polynom: "
					<< cubeSpline(x)<<"; Difference: " <<  diffValue << "\n";
		}
	}
}

template<typename T>
void RunTest::test_HermitePolynomial() {

}

template<typename T>
void RunTest::test_MonotoneCubeHermiteSpline(){

}

template<typename T>
void RunTest::test_CubicPolynomial() {

	const T valueOfPolynom = 2.0;
	CubicPolynomial<T> cpd_constant(1.0,2.0,valueOfPolynom,valueOfPolynom,0.0,0.0);

	//interpolate and check the right value
	if ( std::fabs(cpd_constant(1.5) - valueOfPolynom) > accuracyGoal<T>() ) {
		_allSuccess = false;
		std::cout << "Test of "<< nameOfTypeTrait<T>() << " cubic polynomial failed:"
				"interpolating constant \n";
	}

	//define some test data
	std::vector<T> xValues;
	std::vector<T> dataSet;
	std::vector<T> derivativeData;
	this->get_polynomial_test_data(xValues,dataSet,derivativeData);

	CubicPolynomial<T> cpd_testData(xValues.front(),xValues.back(),
			dataSet.front(),dataSet.back(),
			derivativeData.front(),derivativeData.back());
	for ( size_t i = 0 ; i < xValues.size() ; i++){
		T value = cpd_testData(xValues[i]);
		if ( std::fabs(value - dataSet[i]) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " cubic polynomial failed: "
					"interpolating test data at step " << i << "\n";
			std::cout << "\tData value from polynom: " <<  value << "; Data value test: "
					<< dataSet[i]<<"; Difference: " <<  value - dataSet[i] << "\n";
		}
		T derivative;
		cpd_testData.evaluate_derivative(xValues[i],derivative);
		if ( std::fabs(derivative - derivativeData[i]) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " cubic polynomial failed: "
					"interpolating derivative test data at step " << i << "\n";
			std::cout << "\tDerivative data value from polynom: " <<  derivative << "; Derivative data value test: "
					<< derivativeData[i]<<"; Difference: " <<  derivative - derivativeData[i] << "\n";
		}
	}
};

template<typename T>
void RunTest::get_polynomial_test_data(std::vector<T> & xValues,
		std::vector<T> & polynomialData,
		std::vector<T> & derivativePolynomialData ) const {

	// define the data to pass
	static const T xValuesArray[] = {-0.989218429384234234,-0.2319384234234,
			-0.0012134234876,0.62347289433245,1.3453158523413245};
	static const T dataArray[] = {-5.3891789678248282332,.8209848419213039104,
			0.9997192497522201019,1.9744090986534953532,12.0114397015586715940};
	static const T derivativeData[] = {17.8055531218330275340,1.5949988948483670073,
			0.2327484514620032387,4.9287797622237035407,25.6058553748127175706};

	//	put it into a clean vector
	xValues.clear();
	xValues.insert(xValues.end(),xValuesArray,xValuesArray+sizeof(xValuesArray)/sizeof(T));

	polynomialData.clear();
	polynomialData.insert(polynomialData.end(),dataArray,dataArray+sizeof(dataArray)/sizeof(T));

	derivativePolynomialData.clear();
	derivativePolynomialData.insert(derivativePolynomialData.end(),derivativeData,
			derivativeData+sizeof(derivativeData)/sizeof(T));
}

};/* namespace data_interpolation */
};/* namespace gslpp */

