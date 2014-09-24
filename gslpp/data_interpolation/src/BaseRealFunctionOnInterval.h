/*
 * BasePolynomial.h
 *
 *  Created on: Sep 7, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_BASEPOLYNOMIAL_H_
#define GSLPP_DATA_INTERPOLATION_BASEPOLYNOMIAL_H_

namespace gslpp {
namespace data_interpolation {

/**	A base class for shared code of real valued functions on a real valued interval interval.
 *
 * 	The right (larger) border is not part of the range of definition.
 * 	A polynomial obeys a comparison hirachie. We say one cubic polynomial is smaller than the other
 * 	if the largest value is smaller than smallest value of the respective other polynomial.
 */
template<typename T,class derived>
class BaseRealFunctionOnInterval {
public:
	/**	Constructor setting the range of definition to NaN
	 */
	BaseRealFunctionOnInterval();

	/**	Constructor setting the range of definition to rangeMin and rangeMax.
	 *
	 *	Sets the internal init status to true.
	 *
	 * @param rangeMin
	 * @param rangeMax
	 */
	BaseRealFunctionOnInterval(T rangeMin, T rangeMax);

	/** Evaluate at position x.
	 *
	 * Calls this->evaluate(T x) internally
	 */
	T operator() (T x) const;

	/** Get the infinium of the range of definition
	 *
	 * @return Infinum of the range of definition
	 */
	T min_range() const;

	/** Get the suppremum of the range of definition
	 *
	 * @return Suppremum of the range of definition
	 */
	T max_range() const;

	/**	Get size of the interval.
	 *
	 * @return Size of the interval
	 */
	T interval_length() const;

	/** Checks if x is in the range of definition.
	 *
	 * @param x The position.
	 * @return True if x is in the range of definition, else false.
	 */
	bool x_is_in_range(T x) const;

	/** Check if the polynomial is below the position x.
	 *
	 * @param x The position.
	 * @return True if x is larger or equal than the suppremum of the range of definition
	 */
	bool operator< (T x) const;

	/** Check if the polynomial is above the position x.
	 *
	 * @param x The position.
	 * @return True if x is smaller than the infinum of the range of definition.
	 */
	bool operator> (T x) const;

	/**	Set the infinum and suppremum range of definition to xMin and xMax.
	 *
	 * @param xMin
	 * @param xMax
	 */
	void set_range_of_definition(T xMin, T xMax);

	/**
	 * Overwrites the internal initialization status.
	 *
	 * @param isInit The new status of the object.
	 */
	void set_init_state(bool isInit);

	/**
	 *  Check status of the object.
	 *
	 * @return True if the object is initialized.
	 */
	bool is_init() const;
private:
	bool _isInit;	///If true the object was initialized.
	T _x1; ///Infinum of the range of definition
	T _x2; ///Suppremum of the range of definition
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/BaseRealFunctionOnInterval.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_BASEPOLYNOMIAL_H_ */
