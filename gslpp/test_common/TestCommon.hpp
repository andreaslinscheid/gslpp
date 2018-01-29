/*
 * TestCommon.hpp
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#include "gslpp/test_common/TestCommon.h"
#include "gslpp/auxillary/AccuracyGoal.h"

namespace gslpp {
namespace test_common {

template<typename T>
T TestCommon::accuracyGoal() const {
	return gslpp::auxillary::AccuracyGoal<T>::value;
};

} /* namespace test_common */
} /* namespace gslpp */
