/*
 * Test.hpp
 *
 *  Created on: Sep 5, 2014
 *      Author: alinsch
 */
#include "gslpp/data_interpolation/test/Test.h"
#include "gslpp/data_interpolation/CubicPolynomial.h"
#include "gslpp/data_interpolation/CubeSpline.h"
#include "gslpp/data_interpolation/HermitePolynomial.h"
#include "gslpp/data_interpolation/MonotoneCubeHermiteSpline.h"
#include "gslpp/data_interpolation/BiCubicInterpolation.h"

namespace gslpp{
namespace data_interpolation{

template<typename T>
void RunTest::test_CubeSpline() {

	//Initialize with spline test data that can be solved analytically
	//	(see http://en.wikipedia.org/wiki/Spline_interpolation)
	//This test data leads to two polynomials that are supposed to have the	derivatives
	//	(k0 = -0.6875 , k1 = -0.1250) and (k0 = -0.1250 , k1 = 1.5625 )
	std::vector<T> xValues;
	std::vector<T> dataSet;
	xValues.push_back(-1.0);xValues.push_back(0.0);xValues.push_back(3.0);
	dataSet.push_back( 0.5);dataSet.push_back(0.0);dataSet.push_back(3.0);

	gslpp::data_interpolation::CubeSpline<T,CubicPolynomial<T> > cubeSpline;
	cubeSpline.initialize(xValues,dataSet);

	gslpp::data_interpolation::CubicPolynomial<T> firstPolynom(xValues[0],xValues[1],
			dataSet[0],dataSet[1],-0.6875,-0.1250);
	gslpp::data_interpolation::CubicPolynomial<T> secondPolynom(xValues[1],xValues[2],
			dataSet[1],dataSet[2],-0.1250,1.5625);

	const size_t numEvals = 100;
	for ( size_t i = 0 ; i < numEvals ; i++){

		// sample the interval
		T x = cubeSpline.min_range() + ((cubeSpline.max_range() - cubeSpline.min_range() )/numEvals)*i;

		T diffValue;
		if ( x < xValues[1]) {
			diffValue = firstPolynom(x) - cubeSpline(x);
		} else {
			diffValue = secondPolynom(x) - cubeSpline(x);
		}
		if ( std::fabs(diffValue) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " cubic spline failed: "
					"interpolating test data at step x = " << x << "\n";
			std::cout << "\tData value from spline: " << cubeSpline(x)
					<<"; Data value from polynom: "<< ( x < xValues[1] ? firstPolynom(x) : secondPolynom(x))
					<<"; Difference: " <<  diffValue
					<<"\n\tMismatch at this point with respect to " << ( x < xValues[1] ? "first polynom" : "second polynom")
					<< "\n";
		}
	}
}

template<typename T>
void RunTest::test_HermitePolynomial() {

	//Test for a constant value interpolation
	const T valueOfPolynom = 2.0;
	HermitePolynomial<T> hp_constant(1.0,2.0,valueOfPolynom,valueOfPolynom,0.0,0.0);

	//interpolate and check the right value
	if ( std::fabs(hp_constant(1.5) - valueOfPolynom) > accuracyGoal<T>() ) {
		_allSuccess = false;
		std::cout << "Test of "<< nameOfTypeTrait<T>() << " Hermite polynomial failed:"
				"interpolating constant \n";
	}

	//define some test data
	std::vector<T> xValues;
	std::vector<T> dataSet;
	std::vector<T> derivativeData;
	this->get_polynomial_test_data(xValues,dataSet,derivativeData);

	HermitePolynomial<T> hp_testData(xValues.front(),xValues.back()+std::numeric_limits<T>::epsilon(),
			dataSet.front(),dataSet.back(),
			derivativeData.front(),derivativeData.back());
	for ( size_t i = 0 ; i < xValues.size() ; i++){
		T value = hp_testData(xValues[i]);
		if ( std::fabs(value - dataSet[i]) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " Hermite basis cubic polynomial failed: "
					"interpolating test data at step " << i << "\n";
			std::cout << "\tData value from polynom: " <<  value << "; Data value test: "
					<< dataSet[i]<<"; Difference: " <<  value - dataSet[i] << "\n";
		}
		T derivative;
		hp_testData.evaluate_derivative(xValues[i],derivative);
		if ( std::fabs(derivative - derivativeData[i]) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " Hermite basis cubic polynomial failed: "
					"interpolating derivative test data at step " << i << "\n";
			std::cout << "\tDerivative data value from polynom: " <<  derivative << "; Derivative data value test: "
					<< derivativeData[i]<<"; Difference: " <<  derivative - derivativeData[i] << "\n";
		}
	}
}

template<typename T>
void RunTest::test_MonotoneCubeHermiteSpline(){

	//Initialize with spline test data that can be solved analytically
	//	(see http://en.wikipedia.org/wiki/Monotone_cubic_interpolation)
	//This test data leads to three polynomials that are supposed to have the	derivatives
	//	(k0 = -1 , k1 = 0.0) and (k0 = 0.0 , k1 = 43.0/std::sqrt( 43*43+4.0 ) ) and
	//	(k0 =  43.0/std::sqrt( 43*43+4.0 ) , k1 = 14 )
	std::vector<T> xValues;
	std::vector<T> dataSet;
	xValues.push_back(-2.0);xValues.push_back(-1.0);xValues.push_back(2.0);xValues.push_back( 3.0);
	dataSet.push_back( 1.0);dataSet.push_back( 0.0);dataSet.push_back(1.0);dataSet.push_back(15.0);

	MonotoneCubeHermiteSpline<T,CubicPolynomial<T> > monotoneSpline(xValues,dataSet);

	gslpp::data_interpolation::CubicPolynomial<T> firstPolynom(xValues[0],xValues[1],
			dataSet[0],dataSet[1],-1.0,0.0);
	gslpp::data_interpolation::CubicPolynomial<T> secondPolynom(xValues[1],xValues[2],
			dataSet[1],dataSet[2],0.0, 43.0/std::sqrt( 43*43+4.0 ) );
	gslpp::data_interpolation::CubicPolynomial<T> thirdPolynom(xValues[2],xValues[3],
			dataSet[2],dataSet[3],43.0/std::sqrt( 43*43+4.0 ),14.0);

	const size_t numEvals = 100;
	for ( size_t i = 0 ; i < numEvals ; i++){

		// sample the interval
		T x = monotoneSpline.min_range() + ((monotoneSpline.max_range() - monotoneSpline.min_range() )/numEvals)*i;

		T valueCompare;
		if ( x < xValues[1]) {
			valueCompare = firstPolynom(x);
		} else if ( x < xValues[2] ) {
			valueCompare = secondPolynom(x);
		} else if ( x < xValues[3] ) {
			valueCompare = thirdPolynom(x);
		}
		T diffValue = valueCompare - monotoneSpline(x);

		if ( std::fabs(diffValue) > accuracyGoal<T>() ) {
			_allSuccess = false;
			std::cout << "Test of "<< nameOfTypeTrait<T>() << " monotone Hermite spline failed: "
					"interpolating test data at step x = " << x << "\n";
			std::cout << "\tData value from spline: " << monotoneSpline(x)
					<<"; Data value from polynom: " << valueCompare
					<<"; Difference: " <<  diffValue
					<<"\n\tMismatch at this point with respect to " << ( x < xValues[1] ? "first polynom" :
							( x < xValues[2] ? "second polynom" : "third polynom") )
					<< "\n";
		}
	}
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

	CubicPolynomial<T> cpd_testData(xValues.front(),xValues.back()+std::numeric_limits<T>::epsilon(),
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

template<typename T>
void RunTest::test_CubeSpline2D(){

	//define some test data
	std::vector<T> dataSet;
	std::vector<T> xValues;

	std::vector<T> dataSet1;
	std::vector<T> derivativeData1;
	this->get_polynomial_test_data(xValues,dataSet1,derivativeData1);

	std::vector<T> yValues;
	static const T yValuesArray[] = {-1.0,0.0,1.0};
	yValues.insert(yValues.end(),yValuesArray,yValuesArray+sizeof(yValuesArray)/sizeof(T));

	// preallocate memory and insert identical data along the y axis.
	dataSet.reserve( dataSet1.size()*yValues.size() );
	for ( size_t i = 0 ; i < xValues.size(); i++ ){
		for ( size_t j = 0 ; j < yValues.size(); j++ ){
			dataSet.push_back(dataSet1[i]);
		}
	}

	BiCubicInterpolation<T> cubeSpline2d;
	cubeSpline2d.initialize(xValues,yValues,dataSet);

}

};/* namespace data_interpolation */
};/* namespace gslpp */

