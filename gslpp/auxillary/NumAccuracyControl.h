/*
 * NumAccuracyControl.h
 *
 *  Created on: Oct 20, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_
#define GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_

#include "gslpp/auxillary/has_iterator.h"
#include <cstddef>

namespace gslpp {
namespace auxillary {

template<class derived,typename T, bool THasIterator>
class NumAccuracyControl_impl;

/**
 * 	A class that controls the numerical convergence of an algorithm.
 *
 * 	It allows the comparison with a number of subdivisions of type size_t, local and global
 * 	 absolute and relative error estimation thresholds of type T.
 * 	If T defines an iterator, the elements iterator over must be convertible to a floating point type.
 * 	The behavior is different if T defines an iterator or not. If it does
 * 	 the following applies to all individual elements. Comparisons are performed element wise.
 * 	If a threshold is zero (or all elements are zero) it is not considered.
 * 	The relative comparison is performed by multiplying (all elements of) the function value
 * 	 with (all elements of) the relative error threshold. The result is compared with the estimate for
 * 	 the relative error.
 * 	The criteria of convergence is that all elements are converged.
 */
template<typename T>
class NumAccuracyControl : public NumAccuracyControl_impl<NumAccuracyControl<T>, T, has_iterator<T>::value > {
public:

	/**
	 * 	Constructor using the default settings.
	 *
	 * 	Types that do not define an iterator are initialized with gslpp::auxillary::AccuracyGoal<T>().
	 * 	Other types do not apply default initializations.
	 */
	NumAccuracyControl();

	/**
	 * Check if the local numerical accuracy is sufficient.
	 *
	 * @param localErrEstimate The local error estimate at this point.
	 * @param functionValue The function value at this point.
	 * @return true if localErrEstimate is smaller or equal the local error threshold.
	 */
	bool locally_sufficient(T const& localErrEstimate, T const& functionValue) const;

	/**
	 * Check if the global numerical accuracy is sufficient.
	 *
	 * @param globalErrEstimate The global error estimate.
	 * @param functionValue The global function value to compare with such as the integral.
	 * @return true if globalErrEstimate is smaller or equal the global error threshold.
	 */
	bool global_sufficient(T const& globalErrEstimate, T const& functionValue) const;

	/**
	 * Check if the number of subdivision is below the threshold
	 *
	 * @param numSubdivisons The number of subdivisions.
	 * @return True if \ref numSubdivisons is below the threshold
	 */
	bool sub_divisions_below_max(size_t numSubdivisons) const;

	/**
	 * Performs an internal comparison of first and second.
	 *
	 * For simple types without an iterator, this calls the operator<,
	 * while for types T with an iterator this compares all elements and returns
	 * true if all elements of first are lower than the corresponding element of second.
	 * For dynamic types T, first and second must be of the same size.
	 *
	 * @param first
	 * @param second
	 * @return True if the first is determined to be lower than the second.
	 */
	bool first_lower_than_second(T const& first, T const& second) const;

	/**
	 * Performs an internal comparison of first and second.
	 *
	 * For simple types without an iterator, this calls the operator<,
	 * while for types T with an iterator this compares all elements and returns
	 * true if all elements of first are lower or equal than the corresponding element of second.
	 * For dynamic types T, first and second must be of the same size.
	 *
	 * @param first
	 * @param second
	 * @return True if the first is determined to be lower or equal than the second.
	 */
	bool first_leq_than_second(T const& first, T const& second) const;

	/**
	 * Set the local error threshold to the values given.
	 *
	 * @param thrRel
	 * @param thrAbs
	 */
	void set_local_error_threshold(T const& thrRel, T const& thrAbs);

	/**
	 * Set the subdivision threshold to n.
	 *
	 * @param n The new subdivision threshold
	 */
	void set_max_subdiv(size_t n);

	/**
	 * Set the global error threshold to the values given.
	 *
	 * @param thrRel
	 * @param thrAbs
	 */
	void set_global_error_threshold(T const& thrRel, T const& thrAbs);

	/**
	 * "Overwrite internal subdivision number with n if n is larger.
	 *
	 * @param n The number of subdivisions.
	 */
	void set_num_subdiv(size_t n);

	/**
	 * @return The internal maximal number of subdivision.
	 */
	size_t get_max_num_subdiv() const;

	/**
	 * Tell the object about the absolute error estimate.
	 *
	 * @param estimedAbsErr The absolute error estimate.
	 */
	void set_abs_error_estimate(T const& estimedAbsErr);

	/**
	 * @return The absolute error estimate.
	 */
	T get_abs_error_estimate() const;

	/**
	 * Tell the object about the relative error estimate.
	 *
	 * @param estimedRelErr The relative error estimate.
	 */
	void set_rel_error_estimate(T const& estimedRelErr);

	/**
	 * @return The relative error estimate.
	 */
	T get_rel_error_estimate() const;
private:

	bool _checkAbsLocal;
	bool _checkRelLocal;
	bool _checkAbsGlobal;
	bool _checkRelGlobal;

	T _localRelativeErrorThreshold;
	T _localAbsErrorThreshold;

	T _globalRelativeErrorThreshold;
	T _globalAbsErrorThreshold;

	T _errorEstimateRel;
	T _errorEstimateAbs;

	size_t _subdivisions;
	size_t _maxNumberOfSubdivisions;

};

} /* namespace auxillary */
} /* namespace gslpp */
#include "gslpp/auxillary/src/NumAccuracyControl.hpp"
#endif /* GSLPP_AUXILLARY_NUMACCURACYCONTROL_H_ */
