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
#include <cstddef>
#include "gslpp/data_interpolation/CubicPolynomial.h"

namespace gslpp {
namespace data_interpolation {

/**	A cubic spline that interpolates data.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation]
 */
template<typename T>
class CubeSpline {
public:
	/**	Empty constructor calls just gslpp::data_interpolation::CubeSpline.clear().
	 */
	CubeSpline();

	/**	Constructor that call  gslpp::data_interpolation::CubeSpline.initialize().
	 * @param data
	 * @param mesh
	 */
	CubeSpline(std::vector<T> const& data, std::vector<T> const& mesh);

	/** Evaluate the spine at position x.
	 *
	 * @param x The position.
	 * @param value The value of the polynomial at x.
	 */
	void evaluate(T x, T &value) const;

	/** Evaluate the spine at position x.
	 *
	 * @param x The position.
	 * @param value The value of the polynomial at x.
	 * @param derivative The value of the derivative w.r.t. x of the polynomial at x.
	 */
	void evaluate(T x, T &value, T &derivative) const;

	/** Evaluate the spine at position x.
	 *
	 * @param x The position.
	 * @param value The value of the polynomial at x.
	 * @param derivative The value of the derivative w.r.t. x of the polynomial at x.
	 * @param second_derivative The value of the second derivative w.r.t. x of the polynomial at x.
	 */
	void evaluate(T x, T &value, T &derivative, T &second_derivative) const;

	/**	Erase the content of the spline and set it to the initial state.
	 */
	void clear();

	/**	Initialize the spline.
	 *
	 * @param data A vector with data values i at mesh point i respectively.
	 * @param mesh A vector with position values i at mesh point i respectively.
	 */
	void initialize(std::vector<T> const& data, std::vector<T> const& mesh);

	/** Evaluate the spine at position x.
	 *
	 * Calls gslpp::data_interpolation::CubeSpline<T>::evaluate(T x, T &value) internally.
	 * @param x The position.
	 * @return The value of the polynomial at x.
	 */
	T operator() (T x) const;
private:

	//store the grid values upon initialization
	std::set<T> _gridValuesX;

	//store the polynomials
	std::vector< gslpp::data_interpolation::CubicPolynomial<T> > _polynomials;

	//store interval that was last accessed for increased local evaluation performance
	mutable size_t _lastAccessedPolynomIndex;

	size_t _splineMatrixDim;
	std::vector<T> _splineMatrix;
	std::vector<T> _splineVector;

	size_t find_polynomial_in_range(T x) const;

	void build_spline_matrix( std::vector<T> const& data, std::vector<T> const& mesh );

	void find_derivatives_and_build_polynominals( std::vector<T> const& data, std::vector<T> const& mesh);
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/CubeSpline.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_CUBESPLINE_H_ */
