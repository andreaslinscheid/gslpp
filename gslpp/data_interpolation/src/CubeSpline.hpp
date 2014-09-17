/*
 * CubeSpline.hpp
 *
 *  Created on: Jul 6, 2014
 *      Author: alinsch
 */
#include <cmath>
#include "gslpp/error_handling/Error.h"
#include "gslpp/linear_algebra/LinearAlgebra.h"

namespace gslpp {
namespace data_interpolation {

template<typename T>
CubeSpline<T>::CubeSpline(){
	this->clear();
}

template<typename T>
CubeSpline<T>::CubeSpline(std::vector<T> const& data, std::vector<T> const& mesh){
	this->initialize(data,mesh);
}

template<typename T>
void CubeSpline<T>::clear(){
	_splineMatrix.clear();
	_polynomials.clear();
	_splineVector.clear();
	_lastAccessedPolynomIndex = 0;
}

template<typename T>
void CubeSpline<T>::initialize(std::vector<T> const& data, std::vector<T> const& mesh){
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or ( mesh.size() == 0 ) ) {
		gslpp::error_handling::Error( "Input data for spline generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	this->clear();
	this->build_spline_matrix(data,mesh);
	this->find_derivatives_and_build_polynominals(data,mesh);
}

template<typename T>
void CubeSpline<T>::build_spline_matrix( std::vector<T> const& data, std::vector<T> const& mesh ){
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or ( mesh.size() <= 1 ) ) {
		gslpp::error_handling::Error( "Input data for spline matrix generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif

	//insert the data into the internal mesh
	for ( size_t i = 0 ; i < mesh.size() ; i++ ) {
		_gridValuesX.insert(mesh[i]);
	}

	using std::pow;
	_splineMatrixDim = mesh.size();
	_lastAccessedPolynomIndex = _splineMatrixDim/2;

	//Create the spline matrix
	_splineMatrix.assign(_splineMatrixDim*_splineMatrixDim,0.0);
	for (size_t i=1;i<_splineMatrixDim-1;i++){
		//element i,i-1
		_splineMatrix[i*_splineMatrixDim+i-1] = 1.0/(mesh[i]-mesh[i-1]);
		//element i,i
		_splineMatrix[i*_splineMatrixDim+i] = 2.0*(
				1.0/(mesh[i]-mesh[i-1]) + 1.0/(mesh[i+1]-mesh[i]) );
		//element i,+i
		_splineMatrix[i*_splineMatrixDim+i+1] = 1.0/(mesh[i+1]-mesh[i]);
	}
	//element 0,0
	_splineMatrix[0] = 2.0/(mesh[1]-mesh[0]);
	//element 0,1
	_splineMatrix[1] = 1.0/(mesh[1]-mesh[0]);
	//element N,N-1
	_splineMatrix[(_splineMatrixDim-1)*_splineMatrixDim+_splineMatrixDim-2]
	              = 1.0/(mesh[_splineMatrixDim-1]-mesh[_splineMatrixDim-2]);
	//element N,N
	_splineMatrix[(_splineMatrixDim-1)*_splineMatrixDim+_splineMatrixDim-1]
	              = 2.0/(mesh[_splineMatrixDim-1]-mesh[_splineMatrixDim-2]);

	//create the spline vector
	_splineVector.assign(_splineMatrixDim,0.0);
	_splineVector[0] = 3.0*(data[1]-data[0])/pow(mesh[1]-mesh[0],2);
	for (size_t i=1;i<_splineMatrixDim-1;i++){
		_splineVector[i] = 3.0*((data[i]-data[i-1])/pow(mesh[i]-mesh[i-1],2)
				+(data[i+1]-data[i])/pow(mesh[i+1]-mesh[i],2));
	}
	_splineVector[_splineMatrixDim-1] = 3.0*(data[_splineMatrixDim-1]-data[_splineMatrixDim-2])
			/pow(mesh[_splineMatrixDim-1]-mesh[_splineMatrixDim-2],2);
}

template<typename T>
void CubeSpline<T>::find_derivatives_and_build_polynominals(
		std::vector<T> const& data, std::vector<T> const& mesh){
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or
			( mesh.size()*mesh.size() != _splineMatrix.size() ) or
			( mesh.size() <= 1 ) ) {
		gslpp::error_handling::Error( "Input data for spline matrix generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	gslpp::linear_algebra::LinearAlgebra<T>::invert_square_matrix(_splineMatrix);
	std::vector<T> vectorOfDerivatives(_splineMatrixDim,0.0);
	for (size_t i=0;i<_splineMatrixDim;i++){
		for (size_t j=0;j<_splineMatrixDim;j++)
			vectorOfDerivatives[i] += _splineMatrix[i*_splineMatrixDim + j]*_splineVector[j];
	}
	for (size_t ipol=0;ipol< _splineMatrixDim - 1 ;ipol++){

		// construct a polynomial and insert into the container
		gslpp::data_interpolation::CubicPolynomial<T> cubicPolynom(mesh[ipol],mesh[ipol+1],data[ipol],data[ipol+1],
				vectorOfDerivatives[ipol],vectorOfDerivatives[ipol+1]);
		_polynomials.push_back(cubicPolynom);
	}
}

template<typename T>
size_t CubeSpline<T>::find_polynomial_in_range(T x) const {
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

template<typename T>
void CubeSpline<T>::evaluate(T x, T &value) const {
	_lastAccessedPolynomIndex =  this->find_polynomial_in_range(x);
	 _polynomials[_lastAccessedPolynomIndex].evaluate(x,value);
}

template<typename T>
void CubeSpline<T>::evaluate(T x, T &value, T &derivative) const {
	this->evaluate(x,value);

	//_lastAccessedPolynomIndex is now set to this->find_polynomial_in_range(x)
	//	no need to call it again
	 _polynomials[_lastAccessedPolynomIndex].evaluate_derivative(x,derivative);
}

template<typename T>
void CubeSpline<T>::evaluate(T x, T &value, T &derivative,T &second_derivative) const {
	this->evaluate(x,value,derivative);
	 _polynomials[_lastAccessedPolynomIndex].evaluate_second_derivative(x,second_derivative);
}

template<typename T>
T CubeSpline<T>::operator()(T x) const {
	T value;
	this->evaluate(x,value);
	return value;
};

} /* namespace data_interpolation */
} /* namespace gslpp */
