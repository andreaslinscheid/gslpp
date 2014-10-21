/*
 * Test.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/test/Test.h"
#include "gslpp/integration/Integrator.h"
#include <cmath>
#include <iostream>
#include <assert.h>

namespace gslpp {
namespace integration {

void RunTest::run_test(){
	_allSuccess = true;
	std::cout << "\n\nStarting tests of the integration namespace" <<std::endl;

	test_adaptive_integration<float>();
	test_adaptive_integration<double>();
};

} /* namespace integration */
} /* namespace gslpp */
