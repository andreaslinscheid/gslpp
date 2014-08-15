/*
 * CubeSpline.h
 *
 *  Created on: Jul 6, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_CUBESPLINE_H_
#define GSLPP_DATA_INTERPOLATION_CUBESPLINE_H_

#include <vector>
#include <set>

namespace gslpp {
namespace data_interpolation {
/**	A cubic polynomial.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation].
 * 	The template parameter T is supposed to be float or double.
 * 	The right (larger) border is not part of the range of definition.
 * 	A polynomial obeys a comparison hirachie. We say one cubic polynomial is smaller than the other
 * 	if the largest value is smaller than smallest value of the respective other polynomial.
 */
template<typename T>
class CubicPolynomial {
public:
	/** Constructor that sets the polynomial.
	 *
	 * Internally it calls set_polynomial.
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
	void set_polynomial(T x1,T x2,T y1,T y2,T k1,T k2);

	/**	Evaluate the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @return The value of the polynomial at x.
	 */
	T evaluate(T x) const;

	/** Evaluate the derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @return The value of the derivative of the polynomial at x.
	 */
	T evaluate_derivative(T x) const;

	/** Evaluate the second derivative of the polynomial at x.
	 *
	 *  We allow evaluation outside the range of definition.
	 * @param x The position where to evaluate the polynomial.
	 * @return The value of the second derivative of the polynomial at x.
	 */
	T evaluate_second_derivative(T x) const;

	/** Evaluate at position x.
	 *
	 * Calls gslpp::data_interpolation::CubicPolynomial<T>::evaluate(T x) internally
	 */
	T operator() (T x) const {return this->evaluate(x);};

	/** Get the infinium of the range of definition
	 *
	 * @return Infinum of the range of definition
	 */
	double min() const {return _x1;};

	/** Get the suppremum of the range of definition
	 *
	 * @return Suppremum of the range of definition
	 */
	double max() const {return _x2;};

	/** Checks if x is in the range of definition.
	 *
	 * @param x The position.
	 * @return True if x is in the range of definition, else false.
	 */
	bool x_is_in_range(T x) const { return x>=_x1 and x<_x2;};

	/** Check if the polynomial is below the position x.
	 *
	 * @param x The position.
	 * @return True if x is larger or equal than the suppremum of the range of definition
	 */
	bool operator< (T x) const {return x >= _x2;};

	/** Check if the polynomial is above the position x.
	 *
	 * @param x The position.
	 * @return True if x is smaller than the infinum of the range of definition.
	 */
	bool operator> (T x) const {return x < _x1;};
private:
	T _x1,_x2,_y1,_y2,_a,_b;
};

/**	A cubic spline that interpolates data.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation]
 */
template<typename T>
class CubeSpline {
public:
	CubeSpline<T>::CubeSpline();
private:
	size_t _splineMatrixDim;
	std::vector<T> _splineMatrix;
	std::vector<T> _splineVector;
	mutable size_t _lastAccessedPolynomIndex;
	std::set<gslpp::data_interpolation::CubicPolynomial<T> > _interpolatingPolynomials;
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/CubeSpline.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_CUBESPLINE_H_ */
