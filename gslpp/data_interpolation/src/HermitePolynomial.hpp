/*
 * HermitePolynomial.hpp
 *
 *  Created on: Sep 7, 2014
 *      Author: alinsch
 */

#include "gslpp/error_handling/Error.h"
#include "gslpp/data_interpolation/HermitePolynomial.h"

namespace gslpp{
namespace data_interpolation{

template<typename T>
void HermitePolynomial<T>::initialize(T x1,T x2,T y1,T y2,T k1,T k2){
	_coefficientFAtZero = y1;
	_coefficientFAtOne = y2;

}

template<typename T>
HermitePolynomial<T>::HermitePolynomial(T coeffFAtZero, T coeffFAtOne,
		T coeffDFAtZero,T coeffDFAtOne,
		T xFirst, T xLast)
	:	_coefficientFAtZero( coeffFAtZero ),_coefficientFAtOne(coeffFAtOne ),
		_coefficientdFAtZero( coeffDFAtZero ),_coefficientdFAtOne( coeffDFAtOne ),
		_xFirst( xFirst ),_xLast( xLast ){
	//
	_intervalLength = _xLast - _xFirst;
	//
	if ( _intervalLength < 1e-8 )
		std::Error("Problem in HermitePolynomial: Interval size too small or negative",std::Error::INCORRECT_INIT);
};

template<typename T>
void HermitePolynomial<T>::evaluate(T x, T &fOfX) const {
	T t = (x - _xFirst) /  this->interval_length()l;
	fOfX = _coefficientFAtZero * h00(t) + _coefficientdFAtZero * h10(t) * this->interval_length()
			+ _coefficientFAtOne * h01(t) + _coefficientdFAtOne * h11(t) *  this->interval_length();
}

template<typename T>
void HermitePolynomial<T>::evaluate(T x, T &fOfX, T &dfOfX) const {
	T t = (x - _xFirst) / _intervalLength;
	this->evaluate(x,fOfX);
	dfOfX = (_coefficientFAtZero * dh00(t) + _coefficientdFAtZero * dh10(t) * _intervalLength
			+ _coefficientFAtOne * dh01(t) + _coefficientdFAtOne * dh11(t) * _intervalLength ) / _intervalLength;
}

template<typename T>
void HermitePolynomial<T>::evaluate(T x, T &fOfX, T &dfOfX, T &ddfOfX) const {
	T t = (x - _xFirst) / _intervalLength;
	this->evaluate(x,fOfX,dfOfX);
	ddfOfX = (_coefficientFAtZero * ddh00(t) + _coefficientdFAtZero * ddh10(t) * _intervalLength
			+ _coefficientFAtOne * ddh01(t) + _coefficientdFAtOne * ddh11(t) * _intervalLength )
			 / ( _intervalLength * _intervalLength);
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

}; /* namespace data_interpolation */
}; /* namespace gslpp */
