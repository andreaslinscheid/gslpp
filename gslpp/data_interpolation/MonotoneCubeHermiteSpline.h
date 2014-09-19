/*
 * MonotoneCubeHermiteSpline.h
 *
 *  Created on: Apr 15, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_DATA_INTERPOLATION_MONOTONECUBEHERMITESPLINE_H_
#define GSLPP_DATA_INTERPOLATION_MONOTONECUBEHERMITESPLINE_H_

#include <set>
#include <vector>
#include "gslpp/data_interpolation/HermitePolynomial.h"
#include "gslpp/data_interpolation/src/BaseSpline.h"

namespace gslpp{
namespace data_interpolation {

/**
 * Monotone cubic Hermite interpolation with a continous first derivative using the Fritsch–Carlson method.
 *
 *	The template parameter T is supposed to be float or double.
 *	For details on the algorithm see http://en.wikipedia.org/wiki/Monotone_cubic_interpolation
 */
template<typename T = double, class Polynom = HermitePolynomial<T> >
class MonotoneCubeHermiteSpline : public BaseSpline<MonotoneCubeHermiteSpline<T,Polynom>,T,Polynom> {
public:

	/**
	 * Empty constructor calls just BaseSpline.clear().
	 */
	MonotoneCubeHermiteSpline();

	/**
	 * Constructor that calls MonotoneCubeHermiteSpline.initialize().
	 */
	MonotoneCubeHermiteSpline(std::vector<T> const& mesh, std::vector<T> const& data);

	/**
	 * Initialize the spline.
	 *
	 * @param mesh A vector with position values i at mesh point i respectively.
	 * @param data A vector with data values i at mesh point i respectively.
	 */
	void initialize(std::vector<T> const& mesh,std::vector<T> const& data);

private:

	//compute the derivative values
	void compute_derivatives_no_adjusting(std::vector<T> const& mesh,std::vector<T> const& data,
			std::vector<T> &derivatives, std::vector<T> &secants) const;

	//manipulate derivatives to assure monotonicity
	void ensure_monotonicity_in_derivatives(std::vector<T> &derivatives, std::vector<T> const &secants) const;
};

}; /* namespace data_interpolation */
}; /* namespace gslpp */

#include "gslpp/data_interpolation/src/MonotoneCubeHermiteSpline.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_MONOTONECUBEHERMITESPLINE_H_ */
