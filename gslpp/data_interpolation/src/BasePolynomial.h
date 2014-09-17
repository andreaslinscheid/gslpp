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

/**	A base class for shared code of polynomials.
 *
 */
template<typename T,class derived>
class BasePolynomial {
public:
	/**	Constructor setting the range of definition to NaN
	 */
	BasePolynomial();

	/**	Constructor setting the range of definition to rangeMin and rangeMax.
	 *
	 * @param rangeMin
	 * @param rangeMax
	 */
	BasePolynomial(T rangeMin, T rangeMax);

	/** Evaluate at position x.
	 *
	 * Calls this->evaluate(T x) internally
	 */
	T operator() (T x) const;

	/** Get the infinium of the range of definition
	 *
	 * @return Infinum of the range of definition
	 */
	T min() const;

	/** Get the suppremum of the range of definition
	 *
	 * @return Suppremum of the range of definition
	 */
	T max() const;

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
	 * Sets also _isInit to true.
	 * @param xMin
	 * @param xMax
	 */
	void set_range_of_definition(T xMin, T xMax);

	/** Check status of the object.
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

#include "gslpp/data_interpolation/src/BasePolynomial.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_BASEPOLYNOMIAL_H_ */
