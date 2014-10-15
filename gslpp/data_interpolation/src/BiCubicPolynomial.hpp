/*
 * BiCubicPolynomial.hpp
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/BiCubicPolynomial.h"

namespace gslpp {
namespace data_interpolation {

template<typename T>
BiCubicPolynomial<T>::BiCubicPolynomial(){
	this->clear();
}

template<typename T>
BiCubicPolynomial<T>::BiCubicPolynomial(T f00, T f10, T f01, T f11,
		T dfdx00, T dfdx10, T dfdx01, T dfdx11,
		T dfdy00, T dfdy10, T dfdy01, T dfdy11,
		T d2fdydx00, T d2fdydx10, T d2fdydx01, T d2fdydx11,
		T xMin,T yMin, T xMax, T yMax){
	this->initialize(f00,f10,f01,f11,dfdx00,dfdx10,dfdx01,dfdx11,
			dfdy00,dfdy10,dfdy01,dfdy11,
			d2fdydx00,d2fdydx10,d2fdydx01,d2fdydx11,xMin,yMin,xMax,yMax);
}

template<typename T>
void BiCubicPolynomial<T>::initialize(T f00, T f10, T f01, T f11,
					T dfdx00, T dfdx10, T dfdx01, T dfdx11,
					T dfdy00, T dfdy10, T dfdy01, T dfdy11,
					T d2fdydx00, T d2fdydx10, T d2fdydx01, T d2fdydx11,
					T xMin,T yMin, T xMax, T yMax){

	//conditions on the polynomial as an array
	T functionValues[16] = {f00,f10,f01,f11,dfdx00,dfdx10,dfdx01,dfdx11,
			dfdy00,dfdy10,dfdy01,dfdy11,d2fdydx00,d2fdydx10,d2fdydx01,d2fdydx11};

	T intervalLengthX = xMax-xMin;
	T intervalLengthY = yMax-yMin;

	//scale to the unit interval f(x,y)=fs(x',y') ; x' = (x-xMin)/(xMax-xMin) and y the same
	for (size_t i = 0 ; i < 4 ; i++){
		functionValues[i+4] = functionValues[i+4] /intervalLengthX; // df/dx
		functionValues[i+8] = functionValues[i+8] /intervalLengthY; // df/dy
		functionValues[i+12]= functionValues[i+12]/intervalLengthY/intervalLengthX; // ddf/dydx
	}

	//analytic inverse matrix to determine the coefficients
	int InverseCoefficientMatrix [16][16] ={{1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			   {0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
			   {-3,3,0,0,-2,-1,0,0,0,0,0,0,0,0,0,0},
			   {2,-2,0,0,1,1,0,0,0,0,0,0,0,0,0,0},
			   {0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
			   {0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
			   {0,0,0,0,0,0,0,0,-3,3,0,0,-2,-1,0,0},
			   {0,0,0,0,0,0,0,0,2,-2,0,0,1,1,0,0},
			   {-3,0,3,0,0,0,0,0,-2,0,-1,0,0,0,0,0},
			   {0,0,0,0,-3,0,3,0,0,0,0,0,-2,0,-1,0},
			   {9,-9,-9,9,6,3,-6,-3,6,-6,3,-3,4,2,2,1},
			   {-6,6,6,-6,-3,-3,3,3,-4,4,-2,2,-2,-2,-1,-1},
			   {2,0,-2,0,0,0,0,0,1,0,1,0,0,0,0,0},
			   {0,0,0,0,2,0,-2,0,0,0,0,0,1,0,1,0},
			   {-6,6,6,-6,-4,-2,4,2,-3,3,-3,3,-2,-1,-2,-1},
			   {4,-4,-4,4,2,2,-2,-2,2,-2,2,-2,1,1,1,1}};

	//multiply with conditions to get the coefficients
	for (size_t i = 0 ; i < 16 ; i++){
		_coefficients[i] = 0.0;
		for (size_t j = 0 ; j < 16 ; j++)
			_coefficients[i] += functionValues[j]*InverseCoefficientMatrix[i][j];
	}

	//scale to the original interval
	for (size_t i = 0 ; i < 4 ; i++){
		_coefficients[i*4+1] = _coefficients[i*4+1] /intervalLengthX;
		_coefficients[i*4+2] = _coefficients[i*4+2] /intervalLengthX /intervalLengthX;
		_coefficients[i*4+3] = _coefficients[i*4+3] /intervalLengthX /intervalLengthX /intervalLengthX;
	}
	for (size_t i = 0 ; i < 4 ; i++){
		_coefficients[4+i] = _coefficients[4+i] /intervalLengthY;
		_coefficients[8+i] = _coefficients[8+i] /intervalLengthY /intervalLengthY;
		_coefficients[12+i] = _coefficients[12+i] /intervalLengthY /intervalLengthY /intervalLengthY;
	}

	_minX = xMin;
	_minY = yMin;
	_maxX = xMax;
	_maxY = yMax;

	_isIninitalized = true;
}

template<typename T>
void BiCubicPolynomial<T>::clear(){
	_isIninitalized = false;

	_minX = _minY = _maxX = _maxY = std::numeric_limits<T>::quiet_NaN();
	for (size_t i = 0 ; i < 16 ; i++)
		_coefficients[i]= std::numeric_limits<T>::quiet_NaN();
}

template<typename T>
T BiCubicPolynomial<T>::operator ()(T x, T y) const {
#ifdef DEBUG_BUILD
	if ( not this->is_init() )
		gslpp::error_handling::Error("uninitialized access",
				gslpp::error_handling::Error::ACCESS_WITHOUT_INIT);
	if ( not this->is_in_range(x,y) )
		gslpp::error_handling::Error("evaluation outside the range of definition",
				gslpp::error_handling::Error::OUT_OF_BOUNDS);
#endif
	T dataValue;
	this->evaluate(x,y,dataValue);
	return dataValue;
}

template<typename T>
void BiCubicPolynomial<T>::evaluate(T x, T y, T &dataValue) const {
	x = (x -_minX);
	y = (y -_minY);
	dataValue =	_coefficients[0] + _coefficients[1]*x +
			_coefficients[2]*std::pow(x,2) + _coefficients[3]*std::pow(x,3) +
			   (_coefficients[4] + _coefficients[5]*x +
			_coefficients[6]*std::pow(x,2) + _coefficients[7]*std::pow(x,3))*y +
			   (_coefficients[8] + _coefficients[9]*x +
			_coefficients[10]*std::pow(x,2) + _coefficients[11]*std::pow(x,3))*std::pow(y,2) +
			   (_coefficients[12] + _coefficients[13]*x +
			_coefficients[14]*std::pow(x,2) + _coefficients[15]*std::pow(x,3))*std::pow(y,3);
}

template<typename T>
void BiCubicPolynomial<T>::evaluate_derivative(T x, T y, T &gradX, T &gradY) const {
	T deltaX = _maxX-_minX;
	T deltaY = _maxY-_minY;
	x = (x -_minX) ;
	y = (y -_minY) ;

	//The coefficients are scaled by  1/deltaX with the respective power n of x.
	//	For example _coefficients[1] is scaled by 1/deltaX and _coefficients[10] by (1/deltaX)^2.
	//	Since, here, the power of x is reduced by one we need to multiply by deltaX. Furthermore,
	//	taking derivative to the unit interval  means to cast d/dx -> deltaX*d/dx. Thus we multiply
	//	by deltaX*deltaX. Same with y.
	gradX =( ( _coefficients[ 1] + 2.0*_coefficients[ 2]*x + 3.0*_coefficients[ 3]*x*x)       +
			 ( _coefficients[ 5] + 2.0*_coefficients[ 6]*x + 3.0*_coefficients[ 7]*x*x)*y     +
			 ( _coefficients[ 9] + 2.0*_coefficients[10]*x + 3.0*_coefficients[11]*x*x)*y*y   +
			 ( _coefficients[13] + 2.0*_coefficients[14]*x + 3.0*_coefficients[15]*x*x)*y*y*y ) * deltaX*deltaX  ;

	gradY =( ( _coefficients[ 4] + _coefficients[ 5]*x + _coefficients[ 6]*x*x + _coefficients[ 7]*x*x*x)			+
			 ( _coefficients[ 8] + _coefficients[ 9]*x + _coefficients[10]*x*x + _coefficients[11]*x*x*x)*2.0*y 	+
			 ( _coefficients[12] + _coefficients[13]*x + _coefficients[14]*x*x + _coefficients[15]*x*x*x)*3.0*y*y	)*deltaY*deltaY;
}


template<typename T>
void BiCubicPolynomial<T>::evaluate_second_derivative(T x, T y, T &Jxx, T &Jxy, T &Jyy) const {
	T deltaX = _maxX-_minX;
	T deltaY = _maxY-_minY;
	x = (x -_minX) ;
	y = (y -_minY) ;

	//See comment in BiCubicPolynomial<T>::evaluate_derivative for the occurrence of deltaX*deltaX*deltaX*deltaX.
	Jxx = ( ( 2.0*_coefficients[ 2] + 6.0*_coefficients[ 3]*x)       	+
			( 2.0*_coefficients[ 6] + 6.0*_coefficients[ 7]*x)*y    	+
			( 2.0*_coefficients[10] + 6.0*_coefficients[11]*x)*y*y   	+
			( 2.0*_coefficients[14] + 6.0*_coefficients[15]*x)*y*y*y )	*
			deltaX*deltaX*deltaX*deltaX;

	Jxy = ( ( _coefficients[ 5] + 2.0*_coefficients[ 6]*x + 3.0*_coefficients[ 7]*x*x)     		+
			( _coefficients[ 9] + 2.0*_coefficients[10]*x + 3.0*_coefficients[11]*x*x)*2.0*y   	+
			( _coefficients[13] + 2.0*_coefficients[14]*x + 3.0*_coefficients[15]*x*x)*3.0*y*y )*
			deltaX*deltaX*deltaY*deltaY;

	Jyy = ( ( _coefficients[ 8] + _coefficients[ 9]*x + _coefficients[10]*x*x + _coefficients[11]*x*x*x)*2.0 	 +
			( _coefficients[12] + _coefficients[13]*x + _coefficients[14]*x*x + _coefficients[15]*x*x*x)*6.0*y	)*
			deltaY*deltaY*deltaY*deltaY;
}

template<typename T>
bool BiCubicPolynomial<T>::is_in_range(T pointX, T pointY) const {
	return pointX >= _minX ?
			pointX < _maxX ?
					pointY >= _minY ?
							pointY < _maxY ? true : false
							: false
					:false
			: false;
}

template<typename T>
T BiCubicPolynomial<T>::min_range_x() const {
	return _minX;
};

template<typename T>
T BiCubicPolynomial<T>::min_range_y() const {
	return _minY;
};

template<typename T>
T BiCubicPolynomial<T>::max_range_x() const {
	return _maxX;
};

template<typename T>
T BiCubicPolynomial<T>::max_range_y() const {
	return _maxY;
};

template<typename T>
bool BiCubicPolynomial<T>::is_init() const{
	return _isIninitalized;
}

} /* namespace data_interpolation */
} /* namespace gslpp */
