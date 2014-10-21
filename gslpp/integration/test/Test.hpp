/*
 * Test.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/test/Test.h"
#include "gslpp/integration/Integrator.h"
#include <cmath>
#include <iostream>

namespace gslpp {
namespace integration {

template<typename T>
void RunTest::test_adaptive_integration(){
	std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<":" <<std::endl;

	//Integrate f(x)=x from 0 to 2 with the known result 2 for comparison
	auto identityFunctor = []( T x ){ return x; };
	gslpp::integration::Integrator<T,T, decltype( identityFunctor ) > x_integrator;
	T integral,errEstim;
	x_integrator.integrate(0.0,2.0,identityFunctor,integral,errEstim);
	if ( std::fabs(2.0 - integral) > this->accuracyGoal<T>() ){
		std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<" failed.\n" <<
				" Integral of f(x)=x from 0 to 2 not sufficiently close to 2.0. Difference is: " << 2.0 - integral << "\n"<<
				" Instead we have obtained "<< integral <<" from the adaptive integration routine." <<std::endl;
		_allSuccess = false;
	}

	//Integrate the sine function with known result for comparison
	auto sineFunctor = []( T x ){ return std::sin(x); };
	gslpp::integration::Integrator<T,T, decltype( sineFunctor ) > sinus_integrator;
	sinus_integrator.integrate(0.0,M_PI/2.0,sineFunctor,integral,errEstim);
	if ( std::fabs(1.0 - integral) > this->accuracyGoal<T>() ){
		std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<" failed.\n" <<
				" Integral of sine from 0 to pi not sufficiently close to 1.0. Difference is: " << 1.0 - integral << "\n"<<
				" Instead we have obtained "<< integral <<" from the adaptive integration routine." <<std::endl;
		_allSuccess = false;
	}
}

} /* namespace integration */
} /* namespace gslpp */
