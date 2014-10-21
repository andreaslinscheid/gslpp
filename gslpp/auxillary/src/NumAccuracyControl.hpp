/*
 * NumAccuracyControl.hpp
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#include "gslpp/auxillary/NumAccuracyControl.h"
#include "gslpp/auxillary/AccuracyGoal.h"

namespace gslpp {
namespace auxillary {

template<typename T>
NumAccuracyControl<T>::NumAccuracyControl() : localRelativeErrorThreshold( AccuracyGoal<T>::value ),
												localAbsErrorThreshold(  AccuracyGoal<T>::value ),
												globalRelativeErrorThreshold(  AccuracyGoal<T>::value ),
												globalAbsErrorThreshold(  AccuracyGoal<T>::value ) {
}

template<typename T>
bool NumAccuracyControl<T>::locally_sufficient(T localErrEstimate, T functionValue) const {
	return (localErrEstimate <= localRelativeErrorThreshold * functionValue)
			and ( localErrEstimate <= localAbsErrorThreshold);
}

template<typename T>
bool NumAccuracyControl<T>::global_sufficient(T globalErrEstimate, T functionValue) const {
	return (globalErrEstimate <= globalRelativeErrorThreshold * functionValue)
			and ( globalErrEstimate <= globalAbsErrorThreshold);
}

} /* namespace auxillary */
} /* namespace gslpp */
