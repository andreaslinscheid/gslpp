/*
 * BaseSpline.hpp
 *
 *  Created on: Sep 18, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/src/BaseSpline.h"

namespace gslpp {
namespace data_interpolation {

template<class derived, typename T,class polynom>
BaseSpline<derived,T,polynom>::BaseSpline() {
	BaseSpline<derived,T,polynom>::clear();
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::clear() {
	_polynomials.clear();
	_gridValuesX.clear();
	_lastAccessedPolynomIndex = 0;
	this->set_init_state(false);
}

template<class derived, typename T,class polynom>
BaseSpline<derived,T,polynom>::BaseSpline(std::vector<T> const& mesh,
		std::vector<T> const& data) {
	static_cast<derived*>(this)->initialize(mesh,data);
}

template<class derived, typename T,class polynom>
size_t BaseSpline<derived,T,polynom>::find_polynomial_in_range(T x) const {
	if ( _polynomials[_lastAccessedPolynomIndex].x_is_in_range(x) )
		return _lastAccessedPolynomIndex;
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
	_lastAccessedPolynomIndex = static_cast<size_t>( indexX );
	return _lastAccessedPolynomIndex;
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::evaluate(T x, T &value) const {
	_lastAccessedPolynomIndex =  this->find_polynomial_in_range(x);
	 _polynomials[_lastAccessedPolynomIndex].evaluate(x,value);
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::evaluate(T x, T &value, T &derivative) const {
	this->evaluate(x,value);

	//_lastAccessedPolynomIndex is now set to this->find_polynomial_in_range(x)
	//	no need to call it again
	 _polynomials[_lastAccessedPolynomIndex].evaluate_derivative(x,derivative);
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::evaluate(T x, T &value, T &derivative,T &second_derivative) const {
	this->evaluate(x,value,derivative);
	 _polynomials[_lastAccessedPolynomIndex].evaluate_second_derivative(x,second_derivative);
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::insert_grid(std::vector<T> const& strictlyIncreasingGridX){
#ifdef DEBUG_BUILD
	//check if we have more than 1 point
	if ( strictlyIncreasingGridX.size() <= 1 ) {
		gslpp::error_handling::Error( "Input grid too small" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	_gridValuesX.insert( strictlyIncreasingGridX.begin() , strictlyIncreasingGridX.end() );
#ifdef DEBUG_BUILD
	//check if strictlyIncreasingGridX is sorted
	typename std::set<T>::const_iterator it = _gridValuesX.begin();
	for ( size_t i = 0 ; i < strictlyIncreasingGridX.size() ; ++i){
		if ( *it != strictlyIncreasingGridX[i] ){
			gslpp::error_handling::Error("Input grid not sorted",
					gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
		}
		++it;
	}
#endif
	this->set_range_of_definition(*_gridValuesX.begin(),*_gridValuesX.rbegin());
}

template<class derived, typename T,class polynom>
void BaseSpline<derived,T,polynom>::insert_polynom(polynom const& p) {
	size_t polynomToBeInserted = _polynomials.size();
#ifdef DEBUG_BUILD

	//check if _gridValuesX is set
	if ( _gridValuesX.empty() ){
		gslpp::error_handling::Error("Input grid not set. Set before inserting polynomials",
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}

	//check if the polynomial firts into the present slot
	//	find the index of the grid value that corresponds to the present polynomial p.
	typename std::set<T>::iterator ptrToLBoundX;
	ptrToLBoundX = _gridValuesX.lower_bound( p.min_range() );
	int indexX = std::distance(_gridValuesX.begin(),ptrToLBoundX)-1;
	if ( *ptrToLBoundX == p.min_range() )
		indexX += 1;
	//	This must be polynomToBeInserted
	if ( static_cast<size_t>(indexX) != polynomToBeInserted ) {
		gslpp::error_handling::Error("Input polynom not in the right order or does not match the grid.",
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	_polynomials.push_back(p);
}

template<class derived, typename T,class polynom>
bool BaseSpline<derived,T,polynom>::is_init() const {
	return ((not _gridValuesX.empty())
			and  ( _gridValuesX.size() == (_polynomials.size()+1) )
			and BaseRealFunctionOnInterval<T, BaseSpline<derived,T,polynom> >::is_init());
}

} /* namespace data_interpolation */
} /* namespace gslpp */
