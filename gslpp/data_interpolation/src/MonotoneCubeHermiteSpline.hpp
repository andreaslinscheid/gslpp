/*
 * MonotoneCubeHermiteSpline.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: alinsch
 */

#include "MonotoneCubeHermiteSpline.h"
#include <cmath>
#include <iostream>

namespace dataInterpolation {

template<typename T>
MonotoneCubeHermiteSpline<T>::MonotoneCubeHermiteSpline(
		std::vector<T> const& strictlyIncreasingGridX,
		std::vector<T> const& functionValuesForXGrid) {
	this->init(strictlyIncreasingGridX,functionValuesForXGrid);
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::evaluate(T x, T &fOfX) const {
	_polynomials[ this->find_polynomial_in_range(x) ].evaluate(x,fOfX);
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::evaluate(T x, T &fOfX, T &dfByDxOfX) const {
	_polynomials[ this->find_polynomial_in_range(x) ].evaluate(x,fOfX,dfByDxOfX);
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::evaluate(T x, T &fOfX, T &dfByDxOfX, T &ddfByDxOfX) const {
	_polynomials[ this->find_polynomial_in_range(x) ].evaluate(x,fOfX,dfByDxOfX,ddfByDxOfX);
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::clear() {
	_gridValuesX.clear();
	_polynomials.clear();
}

template<typename T>
bool MonotoneCubeHermiteSpline<T>::out_of_data_range(T x) const {
	return ( (x < (*_gridValuesX.begin() )) or (x >= *(--_gridValuesX.end() )) );
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::data_range(T &min, T &max) const{
	min = *_gridValuesX.begin();
	max = *(--_gridValuesX.end() );
}

template<typename T>
size_t MonotoneCubeHermiteSpline<T>::find_polynomial_in_range(T x) const {
	if ( _polynomials[_xLastAccess].is_in_range(x) )
		return _xLastAccess;
	//
	//the polynomial index is the index of the point that is not smaller than x
	//	we use <algorithm> lower_bound which gives the iterator to the first element that does not compare less.
	//	The result is thus equivalent to the index of the upper point of the polynomial, except equality which is
	//		used in a different convention and we have to account for this
	typename std::set<T>::const_iterator ptrToLBoundX;
	ptrToLBoundX = _gridValuesX.lower_bound(x);
	int indexX = std::distance(_gridValuesX.begin(),ptrToLBoundX)-1;
	if ( *ptrToLBoundX == x)
		indexX += 1;
	//
	_xLastAccess = static_cast<size_t>( indexX );
	return _xLastAccess;
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::init(std::vector<T> const& strictlyIncreasingGridX,
		std::vector<T> const& functionValuesForXGrid) {
	if ( strictlyIncreasingGridX.size() < 2)
		std::Error("Cannot fit spline through less than two points",std::Error::INTERNAL_LOGIC_CHECK_FAILED);
	if ( functionValuesForXGrid.size() != strictlyIncreasingGridX.size() )
		std::Error("#Gridpoints unequal to #Function values",std::Error::INTERNAL_LOGIC_CHECK_FAILED);
	//
	this->clear();
	//
	_gridValuesX.insert( strictlyIncreasingGridX.begin() , strictlyIncreasingGridX.end() );
#ifdef __DEBUG
	//
	//check if strictlyIncreasingGridX is sorted
	typename std::set<T>::const_iterator it = _gridValuesX.begin();
	for ( size_t i = 0 ; i < strictlyIncreasingGridX.size() ; ++i){
		if ( *it != strictlyIncreasingGridX[i] ){
			std::Error("input grid appears to be not sorted",std::Error::INTERNAL_LOGIC_CHECK_FAILED);
		}
		++it;
	}
#endif
	//
	//evaluate derivatives
	std::vector<T> derivatives;
	std::vector<T> secants;
	this->compute_derivatives_no_adjusting( functionValuesForXGrid , derivatives , secants );
	//
	//ensure monotonicity
	ensure_monotonicity_in_derivatives(derivatives,secants);
	//
	//create the polynomials
	_polynomials.reserve( _gridValuesX.size()-1 );
	typename std::set<T>::const_iterator xi = _gridValuesX.begin();
	typename std::set<T>::const_iterator xiPlus1 =  _gridValuesX.begin();
	++xiPlus1;
	for ( size_t i = 0 ; i < _gridValuesX.size()-1 ; ++i ){
		dataInterpolation::HermitePolynomial<T> polynom( functionValuesForXGrid[i], functionValuesForXGrid[i+1],
				derivatives[i], derivatives[i+1],*xi,*xiPlus1);
		_polynomials.push_back( polynom );
		++xi;
		++xiPlus1;
	}
	//
	//set a default last access
	_xLastAccess = _polynomials.size()/2;
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::compute_derivatives_no_adjusting(
		std::vector<T> const& functionValues, std::vector<T> &derivatives, std::vector<T> &secants) const {
	if ( functionValues.size() < 2 )
		std::Error("Cannot fit spline through just one point",std::Error::INTERNAL_LOGIC_CHECK_FAILED);
	if ( functionValues.size() != _gridValuesX.size() )
		std::Error("#Gridpoints unequal to #Function values",std::Error::INTERNAL_LOGIC_CHECK_FAILED);
	//
	derivatives.clear();
	derivatives.reserve(functionValues.size());
	secants.clear();
	secants.reserve(functionValues.size()-1);
	//
	//compute secants for each interval
	typename std::set<T>::const_iterator it = _gridValuesX.begin();
	T interval;
	for ( size_t i = 0 ; i < functionValues.size() - 1 ; ++i){
		interval = *it;
		++it;
		interval = *it - interval;
		secants.push_back( (functionValues[i+1]-functionValues[i])/interval );
	}
	//
	//compute derivatives for each point as the average of the two secants of the adjacent intervals
	//	the first and the last one are set to the secants of that interval
	derivatives.push_back( secants.front() );
	for ( size_t i = 1 ; i < functionValues.size() - 1 ; ++i)
		derivatives.push_back( (secants[i-1] + secants[i])/2.0 );
	derivatives.push_back( secants.back() );
}

template<typename T>
void MonotoneCubeHermiteSpline<T>::ensure_monotonicity_in_derivatives(
		std::vector<T> &derivatives, std::vector<T> const &secants) const {
	T lastBeta = 0;
	for ( size_t i = 0 ; i < _gridValuesX.size() - 1; ++i){
		if ( std::fabs(secants[i]) < 1e-8 ){
			derivatives[i] = derivatives[i+1] = 0;
		} else {
			T alpha = derivatives[i]/secants[i];
			T beta = derivatives[i+1]/secants[i];
			//
			//
			if ( (alpha*alpha + beta*beta) > 9 ) {
				derivatives[i] = 3 *  secants[i] * alpha / std::sqrt(  alpha*alpha + beta*beta );
				derivatives[i+1] = 3 *  secants[i] * beta / std::sqrt(  alpha*alpha + beta*beta );
			}
			//
			if ( ((alpha < 0) or (lastBeta < 0)) )
				derivatives[i] = 0;
			//
			lastBeta = beta;
		}
	}
}

template<typename T>
T MonotoneCubeHermiteSpline<T>::interpolate_with_zero_out_of_range(T x) const {
	return this->out_of_data_range(x) ? 0 : (*this)(x);
}

} /* namespace dataInterpolation */
