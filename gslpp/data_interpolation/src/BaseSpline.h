/*
 * BaseSpline.h
 *
 *  Created on: Sep 18, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_BASESPLINE_H_
#define GSLPP_DATA_INTERPOLATION_BASESPLINE_H_

#include <vector>
#include <set>
#include <cstddef>
#include "gslpp/error_handling/Error.h"
#include "gslpp/data_interpolation/CubicPolynomial.h"

namespace gslpp {
namespace data_interpolation {

/**
 * 	Base spline for shared code.
 */
template<class derived,
		typename T = double,
		class polynom = CubicPolynomial<T> >
class BaseSpline : public BaseRealFunctionOnInterval<T, BaseSpline<derived,T,polynom> > {
public:

	/**	Empty constructor calls just gslpp::data_interpolation::CubeSpline.clear().
	 */
	BaseSpline();

	/**	Constructor that call  gslpp::data_interpolation::CubeSpline.initialize().
	 * @param data
	 * @param mesh
	 */
	BaseSpline(std::vector<T> const& mesh,std::vector<T> const& data);

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

	/**	Check the internal status of the object.
	 *
	 * @return True If the object is correctly initialized.
	 */
	bool is_init() const;
protected:

	/**	Set the internal mesh of x values. Also sets BaseRealFunctionOnInterval.set_range_of_definition .
	 *
	 * @param strictlyIncreasingGridX The mesh of \f$x\f$ values where \f$x_i < x_i+1\f$
	 */
	void insert_grid(std::vector<T> const& strictlyIncreasingGridX);

	/**	Insert the polynoms corresponding to the values in _gridValuesX one after the other.
	 *
	 *	Call after setting the grid.
	 *	On DEBUG_BUILD the method checks if the polynom matches the x mesh and fails if not.
	 *
	 * @param p The polynom used for interpolation in this range.
	 */
	void insert_polynom(polynom const& p);
private:

	//store the grid values upon initialization
	std::set<T> _gridValuesX;

	//store the polynomials
	std::vector<polynom> _polynomials;

	//store interval that was last accessed for increased local evaluation performance
	mutable size_t _lastAccessedPolynomIndex;

	size_t find_polynomial_in_range(T x) const;
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/BaseSpline.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_BASESPLINE_H_ */
