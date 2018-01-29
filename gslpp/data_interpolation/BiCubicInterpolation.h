/*
 * CubeSplineInterpolation2D.h
 *
 *  Created on: Nov 26, 2013
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_CUBESPLINE2D_H_
#define GSLPP_DATA_INTERPOLATION_CUBESPLINE2D_H_

#include <cstdlib>
#include <vector>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iterator>
#include "gslpp/data_interpolation/BiCubicPolynomial.h"

namespace gslpp {
namespace data_interpolation {

/**
 * 	2D Data interpolator using polynomials in between input data.
 *
 * 	Upon initialization a method to generate the conditions for the individual interpolation polynomials can be chosen.
 * 	This may be monotone or not. If monotone is chosen the underlying derivatives are taken from MonontoneCubeHermiteSpline
 * 	while otherwise the object uses the CubeSpline to generate smoother data.
 *
 * 	Internal data layout is with y low to high running faster and x low to high running slower.
 * 	 Thus, upon multiple evaluation, memory locality may benefit from evaluation in this respective order.
 */
template<typename T>
class BiCubicInterpolation {
public:

	/**
	 * Empty constructor calls BiCubicInterpolation.clear()
	 */
	BiCubicInterpolation();

	/**
	 * Constructor that sets the BiCubicInterpolation. Internally calls
	 * BiCubicInterpolation.initialize()
	 */
	BiCubicInterpolation(std::vector<T> const& xGridPointValues,
			std::vector<T> const& yGridPointValues,
			std::vector<T> const& data,
			bool monotone = false);

	/**
	 * 	Erase the content and set the object to the inital state.
	 */
	void clear();

	/**
	 * Evaluate the Interpolator at the Point \f$(x,y)\f$.
	 *
	 * Calls BiCubicInterpolation.evaluate internally.
	 */
	T operator() (T x, T y) const;

	/**
	 * Evaluate the Interpolator at the Point \f$(x,y)\f$.
	 *
	 * @param x The x coordinate.
	 * @param y The y coordinate.
	 * @param dataInterpolation Interpolated data at the Point \f$(x,y)\f$.
	 */
	void evaluate(T x, T y, T &dataInterpolation) const;

	/**
	 * Evaluate the derivative of the Interpolator at the Point \f$(x,y)\f$.
	 *
	 * @param x The x coordinate.
	 * @param y The y coordinate.
	 * @param gradX Interpolated derivative w.r.t. x of the data at the Point \f$(x,y)\f$.
	 * @param gradY Interpolated derivative w.r.t. y of the data at the Point \f$(x,y)\f$.
	 */
	void evaluate_derivative(T x, T y,T &gradX, T &gradY) const;

	/**
	 * Compute the BiCubic interpolator from input data.
	 *
	 * The method is to take the regular grid and create spline along each row and column. Then
	 *	use the continuous derivatives to determine df/dx and df/dy. Finally create another
	 *	set of splines from the df/dx values along y at the grid points where the y derivative
	 *	gives the d^2f/dxdy values.
	 * The spline defined the behavior of the interpolator.
	 * The data must be ordered with x (i) slow and y (j) indices fast running in memory:
	 *		data[i*numPtsY + j]
	 *
	 * @param xGridPointValues The x coordinates of the regular grid values.
	 * @param yGridPointValues The y coordinates of the regular grid values.
	 * @param data Data values in the memory line out \f$[x_i \cdot yGridPointValues.size() + y_j]\f$
	 * @param monotone If false a cubic spline is used to determine the derivatives
	 * 				   which leads to a continuous second derivative. Else the MonotoneHermiteSpline
	 * 				   is used which leads to a strict monontonicity between input data values.
	 */
	void initialize(std::vector<T> const& xGridPointValues,
			std::vector<T> const& yGridPointValues,
			std::vector<T> const& data,
			bool monotone = false);

	/**
	 * @return the infinium of the range of definition in x
	 */
	T min_range_x() const;

	/**
	 * @return the infinium of the range of definition in y
	 */
	T min_range_y() const;

	/**
	 * @return the suppremum of the range of definition in x
	 */
	T max_range_x() const;

	/**
	 * @return the suppremum of the range of definition in y
	 */
	T max_range_y() const;

	/**
	 * Get the range of definition.
	 */
	void data_range(T &xMin, T &xMax, T &yMin, T &yMax) const;
private:

	//store last access to speed up the local multiple evaluation speed
	mutable size_t _xLastAccess, _yLastAccess;

	//
	std::vector<BiCubicPolynomial<T> > _interpolatingPolynomials;
	//
	//number of polynomials in each direction
	size_t _numPolynomsX, _numPolynomsY;
	//
	//the data range that the initialized data covers
	T _minRangeX, _minRangeY, _maxRangeX, _maxRangeY;
	//
	//the grid point x coordinates
	std::set<T> _gridValuesX;
	//
	//the grid point y coordinates
	std::set<T> _gridValuesY;
	//
	//flag true if the object is ready for access
	bool _isInit;
	//
	//find the polynomial index that has argument1,argument 2 in range
	//	call an error it is in the range of no one
	size_t find_polynomial_in_range(T argument1, T argument2) const ;
	//
	//initialize of members to zero
	void set_to_zero();

	template<class Plyonom>
	void initialize_t(std::vector<T> const& xGridPointValues,
				std::vector<T> const& yGridPointValues,
				std::vector<T> const& data);
};



};/* namespace data_interpolation */
};/* namespace gslpp */

#include "gslpp/data_interpolation/src/BiCubicInterpolation.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_CUBESPLINE2D_H_ */
