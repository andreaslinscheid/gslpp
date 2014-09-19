/*
 * BasePolynomial.hpp
 *
 *  Created on: Sep 7, 2014
 *      Author: alinsch
 */
#include <limits>
#include "gslpp/data_interpolation/src/BaseRealFunctionOnInterval.h"
#include "gslpp/error_handling/Error.h"

namespace gslpp {
namespace data_interpolation {

template<typename T, class derived>
T BaseRealFunctionOnInterval<T,derived>::operator() (T x) const{
#ifdef DEBUG_BUILD
	if ( not this->is_init() )
		gslpp::error_handling::Error("uninitialized access",
				gslpp::error_handling::Error::ACCESS_WITHOUT_INIT);
	if ( not this->x_is_in_range(x) )
		gslpp::error_handling::Error("evaluation outside the range of definition",
				gslpp::error_handling::Error::OUT_OF_BOUNDS);
#endif
	T val;
	static_cast<derived const * >(this)->evaluate(x,val);
	return val;
};

template<typename T, class derived>
BaseRealFunctionOnInterval<T,derived>::BaseRealFunctionOnInterval() : _isInit(false) {
	this->set_range_of_definition( std::numeric_limits<T>::quiet_NaN() , std::numeric_limits<T>::quiet_NaN());
};

template<typename T, class derived>
BaseRealFunctionOnInterval<T,derived>::BaseRealFunctionOnInterval(T rangeMin, T rangeMax){
	this->set_range_of_definition(rangeMin,rangeMax);
	this->set_init_state(true);
};

template<typename T, class derived>
void BaseRealFunctionOnInterval<T,derived>::set_range_of_definition(T xMin, T xMax){
	_x1 = xMin;
	_x2 = xMax;
}

template<typename T, class derived>
void BaseRealFunctionOnInterval<T,derived>::set_init_state(bool isInit){
	_isInit = true;
}

template<typename T, class derived>
bool BaseRealFunctionOnInterval<T,derived>::is_init() const {
	return _isInit;
}

template<typename T, class derived>
bool BaseRealFunctionOnInterval<T,derived>::operator< (T x) const {
	return x >= _x2;
};

template<typename T, class derived>
bool BaseRealFunctionOnInterval<T,derived>::operator> (T x) const {
	return x < _x1;
};

template<typename T, class derived>
bool BaseRealFunctionOnInterval<T,derived>::x_is_in_range(T x) const {
	return x>=_x1 and x<_x2;
};

template<typename T, class derived>
T BaseRealFunctionOnInterval<T,derived>::interval_length() const {
	return this->max_range() - this->min_range();
};

template<typename T, class derived>
T BaseRealFunctionOnInterval<T,derived>::min_range() const {
	return _x1;
};

template<typename T, class derived>
T BaseRealFunctionOnInterval<T,derived>::max_range() const {
	return _x2;
};

} /* namespace data_interpolation */
} /* namespace gslpp */
