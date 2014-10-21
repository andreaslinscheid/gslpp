/*
 * NumAccuracyControl.h
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_
#define GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_

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
private:

	T localRelativeErrorThreshold;
	T localAbsErrorThreshold;

	T globalRelativeErrorThreshold;
	T globalAbsErrorThreshold;

};

} /* namespace auxillary */
} /* namespace gslpp */
#include "gslpp/auxillary/src/NumAccuracyControl.hpp"
#endif /* GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_ */
