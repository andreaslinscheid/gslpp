/*
 * CubicPolynomial.hpp
 *
 *  Created on: Aug 28, 2014
 *      Author: alinsch
 */
#include <cmath>
#include "gslpp/data_interpolation/CubicPolynomial.h"

namespace gslpp {
namespace data_interpolation {


template<typename T>
CubicPolynomial<T>::CubicPolynomial(T x1,T x2,T y1,T y2,T k1,T k2) {
	this->initialize(x1,x2,y1,y2,k1,k2);
}

template<typename T>
void CubicPolynomial<T>::initialize(T x1,T x2,T y1,T y2,T k1,T k2){
	this->set_range_of_definition(x1,x2);
	_a = k1*( x2 - x1 ) - ( y2 - y1 );
	_b = -k2*( x2 - x1 ) + ( y2 - y1 );
	_y1 = y1;
	_y2 = y2;
	this->set_init_state(true);
}

template<typename T>
void CubicPolynomial<T>::evaluate(T x, T &value) const {
	T t = ( x - this->min_range() ) / this->interval_length() ;
	value = (1-t)*_y1+t*_y2+t*(1-t)*(_a*(1-t)+_b*t);
}

template<typename T>
void CubicPolynomial<T>::evaluate_derivative(T x, T &value) const {
	T t = ( x - this->min_range() ) / this->interval_length() ;
	value = (_y2-_y1)/this->interval_length()+(1-2*t)*(_a*(1-t)+_b*t)
			/this->interval_length() + t*(1-t)*(_b-_a)/this->interval_length();
}

template<typename T>
void CubicPolynomial<T>::evaluate_second_derivative(T x, T &value) const {
	T t = ( x - this->min_range() ) / this->interval_length();
	value = 2*(_b-2*_a+(_a-_b)*3*t)/std::pow(this->interval_length(),2);
}

template<typename T>
T CubicPolynomial<T>::derivative_at_range_min() const {
	return (_y2-_y1+_a)/this->interval_length();
}

template<typename T>
T CubicPolynomial<T>::derivative_at_range_max() const {
	return (_y2-_y1-_b)/this->interval_length();
}

} /* namespace data_interpolation */
} /* namespace gslpp */
