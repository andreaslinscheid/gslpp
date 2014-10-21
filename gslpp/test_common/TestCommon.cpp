/*
 * TestCommon.hpp
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#include "gslpp/test_common/TestCommon.h"
#include <string>

namespace gslpp {
namespace test_common {

template<>
std::string TestCommon::nameOfTypeTrait<double> () const {
	return "double";
}
template<>
std::string TestCommon::nameOfTypeTrait<float> () const {
	return "float";
}

} /* namespace test_common */
} /* namespace gslpp */
