/*
 * NumAccuracyControl.h
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_
#define GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_

#include <cstddef>

namespace gslpp {
namespace auxillary {

/**
 * 	A class that controls the numerical convergence of an algorithm.
 */
template<typename T>
class NumAccuracyControl {
public:

	/**
	 * 	Constructor using the default settings.
	 */
	NumAccuracyControl();

	/**
	 * Check if the local numerical accuracy is sufficient.
	 *
	 * @param localErrEstimate The local error estimate at this point.
	 * @param functionValue The function value at this point.
	 * @return true if localErrEstimate is smaller or equal the local error threshold.
	 */
	bool locally_sufficient(T localErrEstimate, T functionValue) const;

	/**
	 * Check if the global numerical accuracy is sufficient.
	 *
	 * @param globalErrEstimate The global error estimate.
	 * @param functionValue The global function value to compare with such as the integral.
	 * @return true if globalErrEstimate is smaller or equal the global error threshold.
	 */
	bool global_sufficient(T globalErrEstimate, T functionValue) const;

	bool sub_divisions_below_max(size_t numSubdivisons) const;
private:

	size_t _maxNumberOfSubdivisions;

	T _localRelativeErrorThreshold;
	T _localAbsErrorThreshold;

	T _globalRelativeErrorThreshold;
	T _globalAbsErrorThreshold;

};

} /* namespace auxillary */
} /* namespace gslpp */
#include "gslpp/auxillary/src/NumAccuracyControl.hpp"
#endif /* GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_ */
