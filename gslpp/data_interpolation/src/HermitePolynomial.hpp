/*
 * HermitePolynomial.hpp
 *
 *  Created on: Sep 7, 2014
 *      Author: alinsch
 */
#include "gslpp/data_interpolation/HermitePolynomial.h"

namespace gslpp{
namespace data_interpolation{

template<typename T>
void HermitePolynomial<T>::initialize(T x1,T x2,T y1,T y2,T k1,T k2){
#ifdef DEBUG_BUILD
	//simple check for some input error
	if ( x1 > x2 )
		gslpp::error_handling::Error("Problem in HermitePolynomial:"
				"Interval negative",gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
#endif
	this->set_range_of_definition(x1,x2);

	_coefficientFAtZero = y1;
	_coefficientFAtOne = y2;
	_coefficientdFAtZero = k1;
	_coefficientdFAtOne = k2;

	this->set_init_state(true);

}

template<typename T>
HermitePolynomial<T>::HermitePolynomial()
	:	_coefficientFAtZero( std::numeric_limits<T>::quiet_NaN() ),
	 	_coefficientFAtOne(std::numeric_limits<T>::quiet_NaN() ),
	 	_coefficientdFAtZero( std::numeric_limits<T>::quiet_NaN() ),
	 	_coefficientdFAtOne( std::numeric_limits<T>::quiet_NaN() ){
	this->clear();
}

template<typename T>
HermitePolynomial<T>::HermitePolynomial(T xFirst, T xLast,
		T coeffFAtZero, T coeffFAtOne, T coeffDFAtZero,T coeffDFAtOne ){
	this->initialize(xFirst,xLast,coeffFAtZero,coeffFAtOne,coeffDFAtZero,coeffDFAtOne);
};

template<typename T>
void HermitePolynomial<T>::evaluate(T x, T &value) const {
	T t = (x - this->min_range()) /  this->interval_length();
	value = _coefficientFAtZero * h00(t) + _coefficientdFAtZero * h10(t) * this->interval_length()
			+ _coefficientFAtOne * h01(t) + _coefficientdFAtOne * h11(t) *  this->interval_length();
}

template<typename T>
void HermitePolynomial<T>::evaluate_derivative(T x, T &value) const {
	T t = (x - this->min_range()) / this->interval_length();
	value = (_coefficientFAtZero * dh00(t) + _coefficientdFAtZero * dh10(t) * this->interval_length()
			+ _coefficientFAtOne * dh01(t) + _coefficientdFAtOne * dh11(t) * this->interval_length() ) / this->interval_length();
}

template<typename T>
void HermitePolynomial<T>::evaluate_second_derivative(T x, T &value) const {
	T t = (x - this->min_range()) / this->interval_length();
	value = (_coefficientFAtZero * ddh00(t) + _coefficientdFAtZero * ddh10(t) * this->interval_length()
			+ _coefficientFAtOne * ddh01(t) + _coefficientdFAtOne * ddh11(t) * this->interval_length() )
			 / ( this->interval_length() * this->interval_length());
}

//
//The formulas for the Polynomials
template<typename T>
T HermitePolynomial<T>::h00(T t) const {
	return (1 + 2*t)*(1 - t) * (1 - t);
}
template<typename T>
T HermitePolynomial<T>::h10(T t) const {
	return t*(1 - t)*(1 - t);
}
template<typename T>
T HermitePolynomial<T>::h01(T t) const {
	return t*t*(3 - 2*t);
}
template<typename T>
T HermitePolynomial<T>::h11(T t) const {
	return t*t*(t - 1);
}
//
//The formulas for the derivatives of the Polynomials
template<typename T>
T HermitePolynomial<T>::dh00(T t) const {
	return 6*(t - 1)*t;
}
template<typename T>
T HermitePolynomial<T>::dh10(T t) const {
	return (t - 1)*(3*t - 1);
}
template<typename T>
T HermitePolynomial<T>::dh01(T t) const {
	return -dh00(t);
}
template<typename T>
T HermitePolynomial<T>::dh11(T t) const {
	return t*(3*t - 2);
}
//
//The formulas for the second derivatives of the Polynomials
template<typename T>
T HermitePolynomial<T>::ddh00(T t) const {
	return 12*t - 6;
}
template<typename T>
T HermitePolynomial<T>::ddh10(T t) const {
	return 6*t - 4;
}
template<typename T>
T HermitePolynomial<T>::ddh01(T t) const {
	return - ddh00(t);
}
template<typename T>
T HermitePolynomial<T>::ddh11(T t) const {
	return 6*t - 2;
}


template<typename T>
T HermitePolynomial<T>::derivative_at_range_min() const {
	return _coefficientdFAtZero;
}

template<typename T>
T HermitePolynomial<T>::derivative_at_range_max() const {
	return _coefficientdFAtOne;
}

}; /* namespace data_interpolation */
}; /* namespace gslpp */
