/*
 * Test.cpp
 *
 *  Created on: Oct 15, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/test/Test.h"
#include "gslpp/integration/Integrator.h"
#include "gslpp/auxillary/NumAccuracyControl.h"
#include <cmath>
#include <iostream>

namespace gslpp {
namespace integration {


//a minimal type that return a pair of constant numbers
template<typename T>
class  MinimalObject {
public:
	typedef T value_type;

	MinimalObject() : _val(std::make_pair(0.0,0.0))
	{ };

	MinimalObject(value_type first, value_type second) : _val(std::make_pair(first,second))
	{ };

	MinimalObject & operator* (value_type alpha) {
		_val.first = _val.first*alpha;
		_val.second = _val.second*alpha;
		return *this;
	};

	MinimalObject & operator+ (MinimalObject &rhs) {
		_val.first = _val.first + rhs._val.first;
		_val.second = _val.second + rhs._val.second;
		return *this;
	};

	std::pair<T,T> _val;
};
//a minimal function that returns MinimalObject with x+1 in the first argument and x-1 in the second
template<typename T>
class MinimalFunction {
public:
	MinimalObject<T> operator() (T x) const {
		MinimalObject<T> val(x+static_cast<T>(1.0),x-static_cast<T>(1.0));
		return val;
	};
};


template<typename T>
void RunTest::test_adaptive_integration(){
	std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<":" <<std::endl;

	//Integrate f(x)=x from 0 to 2 with the known result 2 for comparison
	auto identityFunctor = []( T x ){ return x; };
	gslpp::integration::Integrator< decltype( identityFunctor ) > x_integrator;
	T integral;
	gslpp::auxillary::NumAccuracyControl<T> errEstim;
	x_integrator.integrate(0.0,2.0,identityFunctor,integral,errEstim);
	if ( std::fabs(2.0 - integral) > this->accuracyGoal<T>() ){
		std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<" failed.\n" <<
				" Integral of f(x)=x from 0 to 2 not sufficiently close to 2.0. Difference is: " << 2.0 - integral << "\n"<<
				" Instead we have obtained "<< integral <<" from the adaptive integration routine." <<std::endl;
		_allSuccess = false;
	}

	//Integrate the sine function with known result for comparison
	auto sineFunctor = []( T x ){ return std::sin(x); };
	gslpp::integration::Integrator< decltype( sineFunctor ) > sinus_integrator;
	sinus_integrator.integrate(0.0,M_PI/2.0,sineFunctor,integral,errEstim);
	if ( std::fabs(1.0 - integral) > this->accuracyGoal<T>() ){
		std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<" failed.\n" <<
				" Integral of sine from 0 to pi not sufficiently close to 1.0. Difference is: " << 1.0 - integral << "\n"<<
				" Instead we have obtained "<< integral <<" from the adaptive integration routine." <<std::endl;
		_allSuccess = false;
	}

	//Integrate a sharply peaked Lorenzian function
	const T gamma = 0.0001;
	const T x0 = 1.0;
	auto lorenzianFunctor = [&]( T x ){
		return gamma/static_cast<T>(M_PI) / ( static_cast<T>( (x - x0)*(x - x0) ) + gamma*gamma);
	};
	Integrator< decltype( lorenzianFunctor ) > lorenz_integrator;

	lorenz_integrator.integrate(-1000000.0,1000001.0,lorenzianFunctor,integral,errEstim);
	if ( std::fabs(1.0 - integral) > this->accuracyGoal<T>() ){
		std::cout << "\n\tTest of the adaptive integration for type "<< this->nameOfTypeTrait<T>() <<" failed.\n" <<
				" Integral of Lorenz function not sufficiently close to 1.0. Difference is: " << 1.0 - integral << "\n"<<
				" Instead we have obtained "<< integral <<" from the adaptive integration routine." <<std::endl;
		_allSuccess = false;
	}

//	//Integrate a generic minimal object that wraps around a std::pair
//	gslpp::integration::Integrator<MinimalFunction<T> > pair_integrator;
//	MinimalFunction<T> pairWrap;
//	MinimalObject<T> integralPair;
//	gslpp::auxillary::NumAccuracyControl<MinimalObject<T> > errEstimPair;
//	pair_integrator.integrate(-1.0,1.0,pairWrap,integralPair,errEstimPair);
}

} /* namespace integration */
} /* namespace gslpp */
