/*
 * MonotoneCubeHermiteSpline.hpp
 *
 *  Created on: Apr 15, 2014
 *      Author: alinsch
 */

#include "gslpp/data_interpolation/MonotoneCubeHermiteSpline.h"

namespace gslpp {
namespace data_interpolation {

template<typename T,class Polynom>
MonotoneCubeHermiteSpline<T,Polynom>::MonotoneCubeHermiteSpline(){
}

template<typename T,class Polynom>
MonotoneCubeHermiteSpline<T,Polynom>::MonotoneCubeHermiteSpline(
		std::vector<T> const& strictlyIncreasingGridX,
		std::vector<T> const& functionValuesForXGrid)
		: BaseSpline<MonotoneCubeHermiteSpline<T,Polynom>,T,Polynom>(strictlyIncreasingGridX,functionValuesForXGrid){
}

template<typename T,class Polynom>
void MonotoneCubeHermiteSpline<T,Polynom>::initialize(std::vector<T> const& mesh,std::vector<T> const& data) {
#ifdef DEBUG_BUILD
	if ( ( data.size() != mesh.size() ) or ( mesh.size() == 0 ) ) {
		gslpp::error_handling::Error( "Input data for spline generation is rubbish" ,
				gslpp::error_handling::Error::INTERNAL_LOGIC_CHECK_FAILED);
	}
#endif
	this->clear();

	this->insert_grid(mesh);

	std::vector<T> derivatives;
	std::vector<T> secants;
	this->compute_derivatives_no_adjusting( mesh, data , derivatives , secants );

	ensure_monotonicity_in_derivatives(derivatives,secants);

	for ( size_t i = 0 ; i < mesh.size()-1 ; ++i ){
		Polynom polynom( mesh[i],mesh[i+1],data[i], data[i+1],
				derivatives[i], derivatives[i+1]);
		this->insert_polynom(polynom);
	}

	this->set_init_state(true);
}

template<typename T,class Polynom>
void MonotoneCubeHermiteSpline<T,Polynom>::compute_derivatives_no_adjusting(
		std::vector<T> const& mesh,std::vector<T> const& data,
		std::vector<T> &derivatives, std::vector<T> &secants) const {

	derivatives.clear();
	derivatives.reserve(mesh.size());
	secants.clear();
	secants.reserve(mesh.size()-1);

	//compute secants for each interval
	T interval;
	for ( size_t i = 0 ; i < mesh.size() - 1 ; ++i){
		interval = mesh[i+1]-mesh[i];
		secants.push_back( (data[i+1]-data[i])/interval );
	}

	//compute derivatives for each point as the average of the two secants of the adjacent intervals
	//	the first and the last one are set to the secants of that interval
	derivatives.push_back( secants.front() );
	for ( size_t i = 1 ; i < secants.size() ; ++i)
		derivatives.push_back( (secants[i-1] + secants[i])/2.0 );
	derivatives.push_back( secants.back() );
}

template<typename T,class Polynom>
void MonotoneCubeHermiteSpline<T,Polynom>::ensure_monotonicity_in_derivatives(
		std::vector<T> &derivatives, std::vector<T> const &secants) const {
	T lastBeta = 0;
	for ( size_t i = 0 ; i < secants.size() ; ++i){
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


} /* namespace data_interpolation */
} /* namespace gslpp */
