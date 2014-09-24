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
	for ( size_t i = 0 ; i < xGridPointValues.size() ; ++i){
		if ( *it != xGridPointValues[i] )
			gslpp::error_handling::Error("Input grid in x not sorted",
					gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
		++it;
	}
	it = _gridValuesY.begin();
	for ( size_t i = 0 ; i < yGridPointValues.size() ; ++i){
		if ( *it != yGridPointValues[i] )
			gslpp::error_handling::Error("Input grid in y not sorted",
					gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
		++it;
	}
#endif

	//Ensure that we start from a clean state
	this->clear();

	//temporary array with the x derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> xDerivatives(xGridPointValues.size()*yGridPointValues.size(),0.0);

	//generate splines along the x direction
	std::vector<T> sliceOfData(xGridPointValues.size(),0.0);

	//for each y generate an interpolation along x and save the derivatives at the grid points.
	for (size_t iy=0 ; iy < yGridPointValues.size() ; ++iy ){

		for ( size_t ix = 0 ; ix < xGridPointValues.size() ; ++ix)
			sliceOfData[ix] = data[ix*yGridPointValues.size()+iy];

		PolynomialT xPolynomAtIY;
		xPolynomAtIY.initialize(xGridPointValues,sliceOfData);

		//save the derivative values
		for ( size_t ix = 0 ; ix < xGridPointValues.size() ; ++ix){
			T dfYI_dx;
			xPolynomAtIY.evaluate_derivative(xGridPointValues[ix],dfYI_dx);
			xDerivatives[ix*yGridPointValues.size()+iy] = dfYI_dx;
		}
	}

	//temporary array with the y derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> yDerivatives(xGridPointValues.size()*yGridPointValues.size(),0.0);

	//for each x generate an interpolation along y and save the derivatives at the grid points.
	for (size_t ix=0 ; ix < xGridPointValues.size() ; ++ix ){

		for ( size_t iy = 0 ; iy < xGridPointValues.size() ; ++iy)
			sliceOfData[iy] = data[ix*yGridPointValues.size()+iy];

		//use that data[ix*yGridPointValues.size()+0] is the starting point of the y data at this ix
		PolynomialT ySplineAtIX;
		ySplineAtIX.initialize(yGridPointValues, std::vector<T>(
				 &data[ix*yGridPointValues.size()],
				 &data[ix*yGridPointValues.size()]+yGridPointValues.size() ) ); //iterator range constructor

		for ( size_t iy = 0 ; iy < yGridPointValues.size() ; ++iy){
			T dfXI_dy;
			ySplineAtIX.evaluate_derivative(yGridPointValues[iy],dfXI_dy);
			yDerivatives[ix*yGridPointValues.size()+iy] = dfXI_dy;
		}
	}

	//temporary array with the x and y derivatives at the grid values. Has the same ordering as the data.
	std::vector<T> xyDerivatives(xGridPointValues.size()*yGridPointValues.size(),0.0);

	//Interpolate of the derivative of x along the y direction
	for (size_t ix=0 ; ix < xGridPointValues.size() ; ++ix ){

		PolynomialT yPolynomOfDFByDXAtIX;
		yPolynomOfDFByDXAtIX.initialize(yGridPointValues, std::vector<T>(
				&xDerivatives[ix*yGridPointValues.size()],
				&xDerivatives[ix*yGridPointValues.size()]+yGridPointValues.size() ) ); //iterator range constructor

		for ( size_t iy = 0 ; iy < yGridPointValues.size() ; ++iy){
			T d2fdxdy;
			yPolynomOfDFByDXAtIX.evaluate_derivative(yGridPointValues[iy],d2fdxdy);
			xyDerivatives[ix*yGridPointValues.size()+iy] = d2fdxdy;
		}
	}

	//we have yGridPointValues.size()*xGridPointValues.size() points
	//	which makes (yGridPointValues.size()-1)*(xGridPointValues.size()-1) intervals
	_interpolatingPolynomials.assign(
			(yGridPointValues.size()-1)*(xGridPointValues.size()-1),
			BiCubicPolynomial<T>() );

	//
	T cond[16];
	size_t pointIndices[4];
	for ( size_t iIntervalX = 0 ; iIntervalX < (xGridPointValues.size()-1) ; ++iIntervalX ) {
		for ( size_t iIntervalY = 0 ; iIntervalY < (yGridPointValues.size()-1) ; ++iIntervalY ) {

			//set the point indices. We count 1: xmin,ymin 2: xmax,ymin 3: xmin,ymax 4: xmax,ymax
			//	This has to be the same ordering as required by BiCubicPolynomial!
			pointIndices[0] = iIntervalX*yGridPointValues.size() + iIntervalY;
			pointIndices[1] = (iIntervalX+1)*yGridPointValues.size() + iIntervalY;
			pointIndices[2] = iIntervalX*yGridPointValues.size() + (iIntervalY+1);
			pointIndices[3] = (iIntervalX+1)*yGridPointValues.size() + (iIntervalY+1);

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
			_interpolatingPolynomials[iIntervalX*(yGridPointValues.size()-1) + iIntervalY].initialize(
			    		cond[ 0],cond[ 1],cond[ 2],cond[ 3], //data values
			    		cond[ 4],cond[ 5],cond[ 6],cond[ 7], //derivatives along x
			    		cond[ 8],cond[ 9],cond[10],cond[11], //derivatives along y
			    		cond[12],cond[13],cond[14],cond[15], //derivatives along x and y
			    		//set the corner values (x,y)_min and (x,y)_max
			    		xGridPointValues[iIntervalX],yGridPointValues[iIntervalY],
			    		xGridPointValues[iIntervalX+1],yGridPointValues[iIntervalY+1]);
		}
	}

	_numPolynomsX = xGridPointValues.size()-1;
	_numPolynomsY = yGridPointValues.size()-1;
	_minRangeX = *std::min_element(_gridValuesX.begin(),_gridValuesX.end());
	_minRangeY = *std::min_element(_gridValuesY.begin(),_gridValuesY.end());
	_maxRangeX = *std::max_element(_gridValuesX.begin(),_gridValuesX.end());
	_maxRangeY = *std::max_element(_gridValuesY.begin(),_gridValuesY.end());
	_isInit = true;
}

template<typename T>
size_t BiCubicInterpolation<T>::find_polynomial_in_range(T x, T y) const {

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
T BiCubicInterpolation<T>::operator() (T argument1, T argument2) const {
	//
	//first try the last polynomial if it is in range
	if ( _interpolatingPolynomials[_xLastAccess*_numPolynomsY + _yLastAccess].is_in_range(argument1,argument2)){
		T minX=_interpolatingPolynomials[_xLastAccess*_numPolynomsY + _yLastAccess].get_min_x();
		T minY=_interpolatingPolynomials[_xLastAccess*_numPolynomsY + _yLastAccess].get_min_y();
		return _interpolatingPolynomials[_xLastAccess*_numPolynomsY + _yLastAccess](argument1-minX,argument2-minY);
	}
	//
	//see if we are out of the data range in which case we return zero
	if ( (argument1 >= _maxRangeX) or ( argument1 < _minRangeX) or (argument2 >= _maxRangeY) or ( argument2 < _minRangeY)) {
		return 0.0;
	}
	//
	//otherwise find the correct polynomial
	size_t polyIndex = find_polynomial_in_range(argument1,argument2);
	T minX=_interpolatingPolynomials[polyIndex].get_min_x();
	T minY=_interpolatingPolynomials[polyIndex].get_min_y();
	return _interpolatingPolynomials[polyIndex](argument1-minX,argument2-minY);
}

template<typename T>
void BiCubicInterpolation<T>::data_range(T &xMin, T &xMax, T &yMin, T &yMax) const {
	xMin = _minRangeX;
	xMax = _maxRangeX;
	yMin = _minRangeY;
	yMax = _maxRangeY;
}

};/* namespace data_interpolation */
};/* namespace gslpp */
