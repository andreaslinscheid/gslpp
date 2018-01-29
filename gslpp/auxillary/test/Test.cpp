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
#include <array>

namespace gslpp {
namespace auxillary {

//a minimal type that return a pair of constant numbers
template<typename T>
class  MinimalObject : public std::array<T,2> {
public:
	MinimalObject() : std::array<T,2>{{0.0,0.0}}
	{
	};

	MinimalObject(T first, T second) : std::array<T,2>{{first,second}}
	{
	};
};

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

	static_assert(gslpp::auxillary::has_iterator<MinimalObject<float> >::value,
			"has_iterator thinks the MinimalObject<float> type has no iterator, fix that!");
}

} /* namespace auxillary */
} /* namespace gslpp */
