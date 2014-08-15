/*
 * CubeSpline.hpp
 *
 *  Created on: Jul 6, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/CubeSpline.h"

namespace gslpp {
namespace data_interpolation {

template<typename T>
CubicPolynomial<T>::CubicPolynomial(T x1,T x2,T y1,T y2,T k1,T k2) {
	set_polynomial(x1,x2,y1,y2,k1,k2);
}

template<typename T>
void CubicPolynomial<T>::set_polynomial(T x1,T x2,T y1,T y2,T k1,T k2){
	_a = k1*( x2 - x1 ) - ( y2 - y1 );
	_b = -k2*( x2 - x1 ) + ( y2 - y1 );
	_x1 = x1;
	_x2 = x2;
	_y1 = y1;
	_y2 = y2;
}

template<typename T>
T CubicPolynomial<T>::evaluate(T x) const {
	double t = ( x - _x1 ) / ( _x2 - _x1 );
	return (1-t)*_y1+t*_y2+t*(1-t)*(_a*(1-t)+_b*t);
}

template<typename T>
T CubicPolynomial<T>::evaluate_derivative(T x) const {
	double t = ( x - _x1 ) / ( _x2 - _x1 );
	return (_y2-_y1)/(_x2-_x1)+(1-2*t)*(_a*(1-t)+_b*t)/(_x2-_x1)+t*(1-t)*(_b-_a)/(_x2-_x1);
}

template<typename T>
T CubicPolynomial<T>::evaluate_second_derivative(T x) const {
	double t = ( x - _x1 ) / ( _x2 - _x1 );
	return 2*(_b-2*_a+(_a-_b)*3*t)/pow(_x2 - _x1,2);
}
} /* namespace data_interpolation */
} /* namespace gslpp */
