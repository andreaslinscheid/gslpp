/*
 * Test.h
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_INTEGRATION_TEST_H_
#define GSLPP_INTEGRATION_TEST_H_

#include "gslpp/test_common/TestCommon.h"

namespace gslpp {
namespace integration {

class RunTest : public gslpp::test_common::TestCommon {
public:
	void run_test();
private:

	template<typename T>
	void test_adaptive_integration();

	template<typename T>
	void test_non_adaptive_integration();
};

} /* namespace integration */
} /* namespace gslpp */

#include "gslpp/integration/test/Test.hpp"
#endif /* GSLPP_INTEGRATION_TEST_H_ */
