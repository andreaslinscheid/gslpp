/*
 * Test.h
 *
 *  Created on: Jul 6, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_TEST_H_
#define GSLPP_DATA_INTERPOLATION_TEST_H_

#include "gslpp/test_common/TestCommon.h"
#include <cmath>
#include <vector>

namespace gslpp{
namespace data_interpolation{

class RunTest : public gslpp::test_common::TestCommon {
public:
	void run_test();
private:

	template<typename T>
	void test_CubicPolynomial();

	template<typename T>
	void test_CubeSpline();

	template<typename T>
	void test_HermitePolynomial();

	template<typename T>
	void test_MonotoneCubeHermiteSpline();

	template<typename T>
	void test_CubeSpline2D();

	template<typename T>
	void get_polynomial_test_data(std::vector<T> & xValues,
			std::vector<T> & polynomialData,
			std::vector<T> & derivativePolynomialData ) const;
};

}; /* namespace data_interpolation */
}; /* namespace gslpp */

#include "gslpp/data_interpolation/test/Test.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_TEST_H_ */
