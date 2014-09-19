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
#include <cmath>
#include "gslpp/data_interpolation/src/BaseSpline.h"
#include "gslpp/data_interpolation/src/BaseRealFunctionOnInterval.h"
#include "gslpp/error_handling/Error.h"
#include "gslpp/linear_algebra/LinearAlgebra.h"

namespace gslpp {
namespace data_interpolation {

/**
 * A cubic spline that interpolates data.
 *
 * 	The formula is according to Wikipedia [http://en.wikipedia.org/wiki/Spline_interpolation]
 */
template<typename T = double,
		class Polynom = CubicPolynomial<T> >
class CubeSpline : public BaseSpline<CubeSpline<T>,T,Polynom > {
public:

	/**
	 * Empty constructor calls just gslpp::data_interpolation::CubeSpline.clear().
	 */
	CubeSpline();

	/**
	 * Constructor that call  gslpp::data_interpolation::CubeSpline.initialize().
	 *
	 * @param data
	 * @param mesh
	 */
	CubeSpline(std::vector<T> const& data, std::vector<T> const& mesh);

	/**
	 * Initialize the spline.
	 *
	 * @param mesh A vector with position values i at mesh point i respectively.
	 * @param data A vector with data values i at mesh point i respectively.
	 */
	void initialize(std::vector<T> const& mesh,std::vector<T> const& data);

	/**
	 * Removes the content of the object and sets it to its inital state.
	 */
	void clear();
private:

	size_t _splineMatrixDim;
	std::vector<T> _splineMatrix;
	std::vector<T> _splineVector;

	void build_spline_matrix( std::vector<T> const& data, std::vector<T> const& mesh );

	void find_derivatives_and_build_polynominals( std::vector<T> const& data, std::vector<T> const& mesh);
};

} /* namespace data_interpolation */
} /* namespace gslpp */

#include "gslpp/data_interpolation/src/CubeSpline.hpp"
#endif /* GSLPP_DATA_INTERPOLATION_CUBESPLINE_H_ */
