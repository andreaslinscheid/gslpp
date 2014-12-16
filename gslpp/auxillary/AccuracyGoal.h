/*
 * AccuracyGoal.h
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_ACCURACYGOAL_H_
#define GSLPP_AUXILLARY_ACCURACYGOAL_H_
/** \file AccuracyGoal.h
    \brief Define a type trait with numbers for the default accuracy that is required to that type.
*/

#include <complex>
#include <type_traits>

namespace gslpp {
namespace auxillary {

/**
 * 	Specify the number for the default accuracy of the specialized type.
 */
template<typename T>
struct AccuracyGoal {
	static_assert(std::integral_constant<T,false>::value,
			"No specialization for this type T");
};

template<typename T>
struct AccuracyGoal<std::complex<T> > {
	static const std::complex<T> value = std::complex<T>(AccuracyGoal<T>::value,AccuracyGoal<T>::value);
};

template<>
struct AccuracyGoal<double> {
	constexpr static double value = 1e-12;
};

template<>
struct AccuracyGoal<float> {
	constexpr static float value = 1e-5f;
};
} /* namespace auxillary */
} /* namespace gslpp */
#endif /* GSLPP_AUXILLARY_ACCURACYGOAL_H_ */
