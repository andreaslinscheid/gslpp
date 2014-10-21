/*
 * TestCommon.h
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_TEST_COMMON_TESTCOMMON_H_
#define GSLPP_TEST_COMMON_TESTCOMMON_H_

#include "gslpp/auxillary/NumAccuracyControl.h"
#include <string>

namespace gslpp {
namespace test_common {

class TestCommon {
protected:
	bool _allSuccess;
public:
	template<typename T>
	std::string nameOfTypeTrait() const;

	template<typename T>
	T accuracyGoal() const;
};

} /* namespace test_common */
} /* namespace gslpp */
#include "gslpp/test_common/TestCommon.hpp"
#endif /* GSLPP_TEST_COMMON_TESTCOMMON_H_ */
