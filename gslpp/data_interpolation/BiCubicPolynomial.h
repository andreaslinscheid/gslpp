/*
 * BiCubicPolynomial.h
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_BICUBICPOLYNOMIAL_H_
#define GSLPP_DATA_INTERPOLATION_BICUBICPOLYNOMIAL_H_

#include <limits>

namespace gslpp {
namespace data_interpolation {

/**
 * A bi-cubic polynomial on a regular grid.
 *
 */
template<typename T>
class BiCubicPolynomial {
public:
	/**
	 * 	Empty constructor.
	 *
	 * 	Calls BiCubicPolynomial.clear() internally.
	 */
	BiCubicPolynomial();

	/**
	 * 	Initialization constructor.
	 *
	 * 	Calls BiCubicPolynomial.initialize() internally.
	 */
	BiCubicPolynomial(T f00, T f10, T f01, T f11,
			T dfdx00, T dfdx10, T dfdx01, T dfdx11,
			T dfdy00, T dfdy10, T dfdy01, T dfdy11,
			T d2fdydx00, T d2fdydx10, T d2fdydx01, T d2fdydx11,
			T xMin,T yMin, T xMax, T yMax);

	/**
	 * Determine and set the polynomial coefficients.
	 *
	 * Coefficients are matched for the 4 function values, the 4 df/dx values, the 4 df/dy and the 4 d2f/dydx values
	 *	always counting x,y = 0,0 1,0 0,1 1,1 in each category.
	 *	In the last position there are the minimal x,y and then maximal x,y.
	 */
	void initialize(T f00, T f10, T f01, T f11,
			T dfdx00, T dfdx10, T dfdx01, T dfdx11,
			T dfdy00, T dfdy10, T dfdy01, T dfdy11,
			T d2fdydx00, T d2fdydx10, T d2fdydx01, T d2fdydx11,
			T xMin,T yMin, T xMax, T yMax);
	/**
	 * Evaluate the Polynomial at x,y
	 *
	 * This calls BiCubicPolynomial.evaluate() internally
	 */
	T operator() (T x, T y) const;

	/**
	 * Evaluate the Polynomial.
	 *
	 * @param x Position in x.
	 * @param y Position in y.
	 * @param result The Polynomial at this point.
	 */
	void evaluate(T x, T y, T &result) const;

	/**
	 * Evaluate the gradient of the Polynomial.
	 *
	 * @param x Position in x.
	 * @param y Position in y.
	 * @param gradX x-component of the gradient.
	 * @param gradY y-component of the gradient.
	 */
	void evaluate_derivative(T x, T y, T &gradX, T &gradY) const ;

	/**
	 * Evaluate the derivative of the Polynomial.
	 *
	 * @param x Position in x.
	 * @param y Position in y.
	 * @param Jxx xx-component of the Hessian-matrix.
	 * @param Jxy xy-component of the Hessian-matrix.
	 * @param Jyx yx-component of the Hessian-matrix.
	 * @param Jyy yy-component of the Hessian-matrix.
	 */
	void evaluate_second_derivative(T x, T y, T &Jxx, T &Jxy, T &Jyx, T &Jyy) const ;

	/**
	 * Check if a point is in the range of definition
	 *
	 * @param x
	 * @param y
	 * @return true if x,y is in the range of definition
	 */
	bool is_in_range(T x, T y) const;

	/**
	 * @return the infinium of the range of definition in x
	 */
	T get_min_range_x() const;

	/**
	 * @return the infinium of the range of definition in y
	 */
	T get_min_range_y() const;

	/**
	 * @return the suppremum of the range of definition in x
	 */
	T get_max_range_x() const;

	/**
	 * @return the suppremum of the range of definition in y
	 */
	T get_max_range_y() const;

	/**
	 *  Check status of the object.
	 *
	 * @return True if the object is initialized.
	 */
	bool is_init() const;

	/**
	 * Erase the content and set the object to the initial state.
	 */
	void clear();
private:

	///The coefficients of the Polynomial {\sum}_{ij} a(i,j) x ^ i * y ^ j.
	///	Layout i*4+j
	T _coefficients[16];

	///min of the range in x
	T _minX;

	///max of the range in y
	T _maxX;

	///min of the range in y
	T _minY;

	///max of the range in y
	T _maxY;

	///flag if ready for access
	bool _isIninitalized;
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/BiCubicPolynomial.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_BICUBICPOLYNOMIAL_H_ */
