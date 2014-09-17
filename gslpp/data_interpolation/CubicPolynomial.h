/*
 * CubicPolynomial.h
 *
 *  Created on: Aug 28, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/src/BasePolynomial.h"

#ifndef GSLPP_DATA_INTERPOLATION_CUBICPOLYNOMIAL_H_
#define GSLPP_DATA_INTERPOLATION_CUBICPOLYNOMIAL_H_

namespace gslpp {
namespace data_interpolation {

/**	A cubic polynomial for real numbers.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation].
 * 	The template parameter T is supposed to be float or double.
 * 	The right (larger) border is not part of the range of definition.
 * 	A polynomial obeys a comparison hirachie. We say one cubic polynomial is smaller than the other
 * 	if the largest value is smaller than smallest value of the respective other polynomial.
 */
template<typename T>
class CubicPolynomial : public BasePolynomial<T , CubicPolynomial<T> > {
public:
	/** Constructor that sets the polynomial.
	 *
	 * Internally it calls CubicPolynomial.initialize.
	 */
	CubicPolynomial(T x1,T x2,T y1,T y2,T k1,T k2);

	/**	Initialize the Polynomial
	 *
	 * @param x1 The infinum of the range of definition.
	 * @param x2 The suppremum of the range of definition.
	 * @param y1 Data value at x1.
	 * @param y2 Data value at x2.
	 * @param k1 Derivative of data at x1.
	 * @param k2 Derivative of data at x2.
	 */
	void initialize(T x1,T x2,T y1,T y2,T k1,T k2);

	/**	Evaluate the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the polynomial at x.
	 */
	void evaluate(T x, T &value) const;

	/** Evaluate the derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the derivative of the polynomial at x.
	 */
	void evaluate_derivative(T x, T &value) const;

	/** Evaluate the second derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @param value The value of the second derivative of the polynomial at x.
	 */
	void evaluate_second_derivative(T x, T &value) const;

private:
	T _y1; ///data value at _x1
	T _y2; ///data value at _x2
	T _a; ///Derivative parameterization k1*( x2 - x1 ) - ( y2 - y1 )
	T _b; ///Derivative parameterization -k2*( x2 - x1 ) + ( y2 - y1 )
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/CubicPolynomial.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_CUBICPOLYNOMIAL_H_ */
