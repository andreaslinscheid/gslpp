/*
 * MonotoneCubeHermiteSpline.h
 *
 *  Created on: Apr 15, 2014
 *      Author: alinsch
 */

#ifndef MONOTONECUBEHERMITESPLINE_H_
#define MONOTONECUBEHERMITESPLINE_H_

#include "HermitePolynomial.h"
#include <set>
#include <vector>

namespace dataInterpolation {
///
///	Monotone cubic Hermite interpolation using the Fritsch–Carlson method.
///
///	For details on the algorithm see http://en.wikipedia.org/wiki/Monotone_cubic_interpolation .
///
template<typename T>
class MonotoneCubeHermiteSpline {
public:
	//
	//initialize on creation
	MonotoneCubeHermiteSpline(std::vector<T> const& strictlyIncreasingGridX, std::vector<T> const& functionValuesForXGrid);
	//
	//evaluates the interpolating function at x (input)
	T operator() (T x) const {T fOfX; evaluate(x,fOfX); return fOfX; };
	//
	//evaluate the interpolating function at x (input)
	void evaluate (T x, T &fOfX) const;
	//
	//evaluate the interpolating function at x (input) and its derivative
	void evaluate (T x, T &fOfX, T &dfByDxOfX) const;
	//
	//evaluate the interpolating function at x (input) and its first and second derivative
	void evaluate (T x, T &fOfX, T &dfByDxOfX, T &ddfByDxOfX) const;
	//
	//return true if x is outside the data range
	bool out_of_data_range(T x) const;
	//
	//return the data range. min (max) is set to the minimum(maximum) of the stored grid.
	void data_range(T &min, T &max) const;
	//
	//special interpolation function that returns zero if out_of_data_range(x) == true
	T interpolate_with_zero_out_of_range(T x) const;
private:
	//
	//store the grid values upon initialization
	std::set<T> _gridValuesX;
	//
	//store the Hermite polynomials
	std::vector<dataInterpolation::HermitePolynomial<T> > _polynomials;
	//
	//store interval that was last accessed for increased local evaluation performance
	mutable size_t _xLastAccess;
	//
	//compute the spline
	void init(std::vector<T> const& strictlyIncreasingGridX, std::vector<T> const& functionValuesForXGrid);
	//
	//compute the derivative values
	void compute_derivatives_no_adjusting(std::vector<T> const& functionValues, std::vector<T> &derivatives, std::vector<T> &secants) const;
	//
	//manipulate derivatives to assure monotonicity
	void ensure_monotonicity_in_derivatives(std::vector<T> &derivatives, std::vector<T> const &secants) const;
	//
	//clear the content of this object
	void clear();
	//
	//determine the index of the polynomial where x >= _xFirst and x < _xLast
	size_t find_polynomial_in_range(T x) const;
};

} /* namespace dataInterpolation */

#include "MonotoneCubeHermiteSpline.hpp"
#endif /* MONOTONECUBEHERMITESPLINE_H_ */
