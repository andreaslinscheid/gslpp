/*
 * AccuracyGoal.h
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_ACCURACYGOAL_H_
#define GSLPP_AUXILLARY_ACCURACYGOAL_H_

#include <complex>
#include <type_traits>

namespace gslpp {
namespace auxillary {

template<typename T>
struct AccuracyGoal {};

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
