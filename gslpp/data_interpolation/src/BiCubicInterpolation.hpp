/*
 * CubeSplineInterpolation2D.hpp
 *
 *  Created on: Nov 26, 2013
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/BiCubicInterpolation.h"

namespace gslpp {
namespace data_interpolation {

template<typename T>
BiCubicInterpolation<T>::BiCubicInterpolation() {
	this->clear();
}

template<typename T>
BiCubicInterpolation<T>::BiCubicInterpolation(std::vector<T> const& xGridPointValues,
		std::vector<T> const& yGridPointValues,
		std::vector<T> const& data,
		bool monotone) {
	this->initialize(xGridPointValues,yGridPointValues,data,monotone);
}


template<typename T>
void BiCubicInterpolation<T>::clear() {
	_interpolatingPolynomials.clear();
	_gridValuesX.clear();
	_gridValuesY.clear();
	_isInit = false;
	_xLastAccess = 0;
	_yLastAccess = 0;
	_numPolynomsX = 0;
	_numPolynomsY = 0;
}

template<typename T>
void BiCubicInterpolation<T>::initialize(std::vector<T> const& xGridPointValues,
		std::vector<T> const& yGridPointValues,
		std::vector<T> const& data,
		bool monotone){
	if ( monotone ) {
		this->initialize_t<MonotoneCubeHermiteSpline<T> >(xGridPointValues,yGridPointValues,data);
	} else {
		this->initialize_t<CubeSpline<T> >(xGridPointValues,yGridPointValues,data);
	}
}

template<typename T>
template<class PolynomialT>
void BiCubicInterpolation<T>::initialize_t(std::vector<T> const& xGridPointValues,
		std::vector<T> const& yGridPointValues,
		std::vector<T> const& data) {

	const size_t numPtsX = xGridPointValues.size();
	const size_t numPtsY = yGridPointValues.size();

	//Ensure that we start from a clean state
	this->clear();

	//some input checks
	if ( (xGridPointValues.size() < 2) or (yGridPointValues.size() < 2) )
		gslpp::error_handling::Error("initialize with at least 2 points in each direction",
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	if( xGridPointValues.size() * yGridPointValues.size() != data.size() )
		gslpp::error_handling::Error("data size does not match the number of points provided",
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);

	//insert the grids and in case of a debug build check if they were sorted.
	_gridValuesX.insert(xGridPointValues.begin(),xGridPointValues.end());
	_gridValuesY.insert(yGridPointValues.begin(),yGridPointValues.end());
#ifdef DEBUG_BUILD
	typename std::set<T>::const_iterator it = _gridValuesX.begin();
	for ( size_t i = 0 ; i <numPtsX ; ++i){
		if ( *it != xGridPointValues[i] )
			gslpp::error_handling::Error("Input grid in x not sorted",
					gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
		++it;
	}
	it = _gridValuesY.begin();
	for ( size_t i = 0 ; i < numPtsY ; ++i){
		if ( *it != yGridPointValues[i] )
			gslpp::error_handling::Error("Input grid in y not sorted",
					gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
		++it;
	}
#endif

	//temporary array with the x derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> xDerivatives(xGridPointValues.size()*numPtsY,0.0);

	//generate splines along the x direction
	std::vector<T> sliceOfData(xGridPointValues.size(),0.0);

	//for each y generate an interpolation along x and save the derivatives at the grid points.
	for (size_t iy=0 ; iy < numPtsY ; iy++ ){

		for ( size_t ix = 0 ; ix <numPtsX ; ix++)
			sliceOfData[ix] = data[ix*numPtsY+iy];

		PolynomialT xPolynomAtIY;
		xPolynomAtIY.initialize(xGridPointValues,sliceOfData);

		std::vector<T> tmp = xPolynomAtIY.deriviatives_at_underlying_grid_points();

		//save the derivative values
		for ( size_t ix = 0 ; ix <numPtsX ; ix++){
			xDerivatives[ix*numPtsY+iy] = tmp[ix];
		}
	}
	sliceOfData.clear();

	//temporary array with the y derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> yDerivatives(xGridPointValues.size()*numPtsY,0.0);

	sliceOfData.assign(numPtsY,0.0);

	//for each x generate an interpolation along y and save the derivatives at the grid points.
	for (size_t ix=0 ; ix <numPtsX ; ix++ ){

		for ( size_t iy = 0 ; iy <numPtsY ; iy++)
			sliceOfData[iy] = data[ix*numPtsY+iy];

		//use that data[ix*numPtsY+0] is the starting point of the y data at this ix
		PolynomialT ySplineAtIX;
		ySplineAtIX.initialize(yGridPointValues,sliceOfData);

		std::vector<T> tmp = ySplineAtIX.deriviatives_at_underlying_grid_points();

		for ( size_t iy = 0 ; iy < numPtsY ; ++iy){
			yDerivatives[ix*numPtsY+iy] = tmp[iy];
		}
	}

	//temporary array with the x and y derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> xyDerivatives(xGridPointValues.size()*numPtsY,0.0);

	//Interpolate of the derivative of x along the y direction
	for (size_t ix=0 ; ix <numPtsX ; ++ix ){

		//overwrites sliceOfData which has still the correct size
		for ( size_t iy = 0 ; iy <numPtsY ; ++iy)
			sliceOfData[iy] = xDerivatives[ix*numPtsY+iy];

		PolynomialT yPolynomOfDFByDXAtIX;
		yPolynomOfDFByDXAtIX.initialize(yGridPointValues,sliceOfData);

		std::vector<T> tmp = yPolynomOfDFByDXAtIX.deriviatives_at_underlying_grid_points();

		for ( size_t iy = 0 ; iy < numPtsY ; ++iy){
			xyDerivatives[ix*numPtsY+iy] = tmp[iy];
		}
	}

	//we have numPtsY*xGridPointValues.size() points
	//	which makes (numPtsY-1)*(xGridPointValues.size()-1) intervals
	const size_t numIntervalsX =numPtsX-1;
	const size_t numIntervalsY = numPtsY-1;
	_interpolatingPolynomials.assign(numIntervalsX*numIntervalsY, BiCubicPolynomial<T>() );

	T cond[16];
	size_t pointIndices[4];
	for ( size_t iIntervalX = 0 ; iIntervalX < numIntervalsX ; iIntervalX++ ) {
		for ( size_t iIntervalY = 0 ; iIntervalY < numIntervalsY ; iIntervalY++ ) {

			//set the point indices. We count 1: xmin,ymin 2: xmax,ymin 3: xmin,ymax 4: xmax,ymax
			//	This has to be the same ordering as required by BiCubicPolynomial!
			pointIndices[0] = iIntervalX*numPtsY + iIntervalY;
			pointIndices[1] = (iIntervalX+1)*numPtsY + iIntervalY;
			pointIndices[2] = iIntervalX*numPtsY + (iIntervalY+1);
			pointIndices[3] = (iIntervalX+1)*numPtsY + (iIntervalY+1);

			//set the polynomial condition values
			for (size_t indexPT = 0 ; indexPT < 4 ; ++indexPT)
				cond[indexPT] = data[pointIndices[indexPT]];
			for (size_t indexPT = 0 ; indexPT < 4 ; ++indexPT)
				cond[indexPT+4] = xDerivatives[pointIndices[indexPT]];
			for (size_t indexPT = 0 ; indexPT < 4 ; ++indexPT)
				cond[indexPT+8] = yDerivatives[pointIndices[indexPT]];
			for (size_t indexPT = 0 ; indexPT < 4 ; ++indexPT)
				cond[indexPT+12] = xyDerivatives[pointIndices[indexPT]];

			//compute the polynomial
			_interpolatingPolynomials[iIntervalX*(numPtsY-1) + iIntervalY].initialize(
			    		cond[ 0],cond[ 1],cond[ 2],cond[ 3], //data values
			    		cond[ 4],cond[ 5],cond[ 6],cond[ 7], //derivatives along x
			    		cond[ 8],cond[ 9],cond[10],cond[11], //derivatives along y
			    		cond[12],cond[13],cond[14],cond[15], //derivatives along x and y
			    		//set the corner values (x,y)_min and (x,y)_max
			    		xGridPointValues[iIntervalX],yGridPointValues[iIntervalY],
			    		xGridPointValues[iIntervalX+1],yGridPointValues[iIntervalY+1]);
		}
	}

	_numPolynomsX =numPtsX-1;
	_numPolynomsY = numPtsY-1;
	_minRangeX = *std::min_element(_gridValuesX.begin(),_gridValuesX.end());
	_minRangeY = *std::min_element(_gridValuesY.begin(),_gridValuesY.end());
	_maxRangeX = *std::max_element(_gridValuesX.begin(),_gridValuesX.end());
	_maxRangeY = *std::max_element(_gridValuesY.begin(),_gridValuesY.end());
	_isInit = true;
}

template<typename T>
size_t BiCubicInterpolation<T>::find_polynomial_in_range(T x, T y) const {

	//first try the last polynomial if it is in range
	if ( _interpolatingPolynomials[_xLastAccess*_numPolynomsY + _yLastAccess].is_in_range(x,y)){
		return _xLastAccess*_numPolynomsY + _yLastAccess;
	}

	//the polynomial index in x and y is the index of the point that is not smaller than x or y
	//	we use <algorithm> lower_bound which gives the iterator to the first element that does not compare less.
	//	The result is thus equivalent to the index of the upper point of the polynomial, except equality which is
	//		used in a different convention and we have to account for this
	typename std::set<T>::const_iterator ptrToLBoundX;
	ptrToLBoundX = _gridValuesX.lower_bound(x);
	int indexX = std::distance(_gridValuesX.begin(),ptrToLBoundX)-1;
	if ( *ptrToLBoundX == x)
		indexX += 1;
	typename std::set<T>::const_iterator ptrToLBoundY;
	ptrToLBoundY = _gridValuesY.lower_bound(y);
	int indexY = std::distance(_gridValuesY.begin(),ptrToLBoundY)-1;
	if ( *ptrToLBoundY == y)
		indexY += 1;
	return static_cast<size_t>(indexX)*_numPolynomsY+static_cast<size_t>(indexY);
}

template<typename T>
T BiCubicInterpolation<T>::operator() (T x, T y) const {
	T result;
	this->evaluate(x,y,result);
	return result;
}

template<typename T>
void BiCubicInterpolation<T>::evaluate(T x, T y,T &dataInterpolation) const {
	return _interpolatingPolynomials[find_polynomial_in_range(x,y)].evaluate(x,y,dataInterpolation);
}

template<typename T>
void BiCubicInterpolation<T>::evaluate_derivative(T x, T y,T &gradX, T &gradY) const {
	return _interpolatingPolynomials[find_polynomial_in_range(x,y)].evaluate_derivative(x,y,gradX,gradY);
}

template<typename T>
void BiCubicInterpolation<T>::data_range(T &xMin, T &xMax, T &yMin, T &yMax) const {
	xMin = _minRangeX;
	xMax = _maxRangeX;
	yMin = _minRangeY;
	yMax = _maxRangeY;
}


template<typename T>
T BiCubicInterpolation<T>::min_range_x() const {
	return _minRangeX;
};

template<typename T>
T BiCubicInterpolation<T>::min_range_y() const {
	return _minRangeY;
};

template<typename T>
T BiCubicInterpolation<T>::max_range_x() const {
	return _maxRangeX;
};

template<typename T>
T BiCubicInterpolation<T>::max_range_y() const {
	return _maxRangeY;
};

};/* namespace data_interpolation */
};/* namespace gslpp */
