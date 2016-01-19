/*
 * Test.cpp
 *
 *  Created on: Jan 18, 2016
 *      Author: alinsch
 */

#include "Test.h"
#include <set>
#include <complex>
#include <vector>
#include <gslpp/auxillary/has_iterator.h>

namespace gslpp {
namespace auxillary {

void RunTest::run_test()
{
	static_assert(! gslpp::auxillary::has_iterator<double>::value,
			"has_iterator thinks the double type has an iterator, fix that!");

	static_assert(! gslpp::auxillary::has_iterator<std::complex<double> >::value,
			"has_iterator thinks the std::complex<double> type has an iterator, fix that!");

	static_assert(gslpp::auxillary::has_iterator<std::vector<double> >::value,
			"has_iterator thinks the std::vector<double> type has no iterator, fix that!");

	static_assert(gslpp::auxillary::has_iterator<std::set<int>>::value,
			"has_iterator thinks the std::set<double> type has no iterator, fix that!");
}

} /* namespace auxillary */
} /* namespace gslpp */
