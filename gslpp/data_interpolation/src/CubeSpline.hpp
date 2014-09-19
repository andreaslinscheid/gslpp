/*
 * CubeSpline.hpp
 *
 *  Created on: Jul 6, 2014
 *      Author: alinsch
 */
#include "gslpp/data_interpolation/CubeSpline.h"

namespace gslpp {
namespace data_interpolation {

template<typename T, class Polynom>
CubeSpline<T,Polynom>::CubeSpline() : BaseSpline<CubeSpline<T>,T,Polynom>(){
}

template<typename T, class Polynom>
CubeSpline<T,Polynom>::CubeSpline(std::vector<T> const& data, std::vector<T> const& mesh)
	: BaseSpline<CubeSpline<T>,T,Polynom>() {
}

template<typename T, class Polynom>
void CubeSpline<T,Polynom>::clear(){
	_splineMatrix.clear();
	_splineVector.clear();
	BaseSpline<CubeSpline<T>,T,Polynom>::clear();
}

template<typename T, class Polynom>
void CubeSpline<T,Polynom>::initialize(std::vector<T> const& mesh, std::vector<T> const& data){
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or ( mesh.size() == 0 ) ) {
		gslpp::error_handling::Error( "Input data for spline generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	this->clear();
	this->set_range_of_definition(mesh.front(),mesh.back());
	this->build_spline_matrix(data,mesh);
	this->find_derivatives_and_build_polynominals(data,mesh);
	this->set_init_state(true);
}

template<typename T, class Polynom>
void CubeSpline<T,Polynom>::build_spline_matrix( std::vector<T> const& data, std::vector<T> const& mesh ){
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or ( mesh.size() <= 1 ) ) {
		gslpp::error_handling::Error( "Input data for spline matrix generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif

	//insert the data into the internal mesh
	this->insert_grid(mesh);

	using std::pow;
	_splineMatrixDim = mesh.size();

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

template<typename T, class Polynom>
void CubeSpline<T,Polynom>::find_derivatives_and_build_polynominals(
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
		Polynom p(mesh[ipol],mesh[ipol+1],data[ipol],data[ipol+1],
				vectorOfDerivatives[ipol],vectorOfDerivatives[ipol+1]);
		this->insert_polynom(p);
	}
}

} /* namespace data_interpolation */
} /* namespace gslpp */
