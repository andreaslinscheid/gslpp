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
NumAccuracyControl<T>::NumAccuracyControl(){
	 _maxNumberOfSubdivisions =1000;
	_localRelativeErrorThreshold = ( AccuracyGoal<T>::value );
	_localAbsErrorThreshold =(  AccuracyGoal<T>::value );
	_globalRelativeErrorThreshold =(   AccuracyGoal<T>::value );
	_globalAbsErrorThreshold =(  AccuracyGoal<T>::value );
}

template<typename T>
bool NumAccuracyControl<T>::locally_sufficient(T localErrEstimate, T functionValue) const {
	return (localErrEstimate <= _localRelativeErrorThreshold * functionValue)
			and ( localErrEstimate <= _localAbsErrorThreshold);
}

template<typename T>
bool NumAccuracyControl<T>::global_sufficient(T globalErrEstimate, T functionValue) const {
	return (globalErrEstimate <= _globalRelativeErrorThreshold * functionValue)
			and ( globalErrEstimate <= _globalAbsErrorThreshold);
}

template<typename T>
bool NumAccuracyControl<T>::sub_divisions_below_max(size_t numOfSubdivisions) const {
	return _maxNumberOfSubdivisions >= numOfSubdivisions;
}

} /* namespace auxillary */
} /* namespace gslpp */
