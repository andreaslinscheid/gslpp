/*
 * HermitePolynomial.h
 *
 *  Created on: Sep 7, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/src/BaseRealFunctionOnInterval.h"
#include "gslpp/error_handling/Error.h"

#ifndef GSLPP_DATA_INTERPOLATION_HERMITEPOLYNOMIAL_H_
#define GSLPP_DATA_INTERPOLATION_HERMITEPOLYNOMIAL_H_

namespace gslpp{
namespace data_interpolation{

/**	A cubic Hermite polynomial for real numbers.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Cubic_Hermite_spline].
 * 	The template parameter T is supposed to be float or double.
 */
template<typename T>
class HermitePolynomial : public BaseRealFunctionOnInterval<T, HermitePolynomial<T> > {
public:
	/** Empty constructor. Sets the coefficients of the polynomial to NaN.
	 */
	HermitePolynomial();

	/** Constructor that sets the polynomial.
	 *
	 * Internally calls HermitePolynomial.initialize.
	 */
	HermitePolynomial(T x1, T x2,T y1,T y2, T k1, T k2);

	/**	Initialize the Polynomial
	 *
	 * @param x1 The infinum of the range of definition.
	 * @param x2 The suppremum of the range of definition.
	 * @param y1 Data value at x1.
	 * @param y2 Data value at x2.
	 * @param k1 Derivative of data at x1.
	 * @param k2 Derivative of data at x2.
	 */
	void initialize(T x1, T x2,T y1, T y2 ,T k1 ,T k2);

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

	/**
	 * @return The derivative of the data at range infinium.
	 */
	T derivative_at_range_min() const;

	/**
	 * @return The derivative of the data at range suppremum.
	 */
	T derivative_at_range_max() const;
private:

	/** Store the coefficients of the Hermite basis
	 */
	T _coefficientFAtZero; 	///Coefficient of h00
	T _coefficientFAtOne;	///Coefficient of h01
	T _coefficientdFAtZero; ///Coefficient of h10
	T _coefficientdFAtOne;	///Coefficient of h11

	/** Hermite basis functions from http://en.wikipedia.org/wiki/Cubic_Hermite_spline
	 */
	T h00(T t) const;
	T h10(T t) const;
	T h01(T t) const;
	T h11(T t) const;

	/** derivatives of the Hermite basis functions
	 */
	T dh00(T t) const;
	T dh10(T t) const;
	T dh01(T t) const;
	T dh11(T t) const;

	/** second derivatives of the Hermite basis functions
	 */
	T ddh00(T t) const;
	T ddh10(T t) const;
	T ddh01(T t) const;
	T ddh11(T t) const;

};

}; /* namespace data_interpolation */
}; /* namespace gslpp */

#include "gslpp/data_interpolation/src/HermitePolynomial.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_HERMITEPOLYNOMIAL_H_ */
