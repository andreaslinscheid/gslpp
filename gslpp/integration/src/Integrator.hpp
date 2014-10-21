/*
 * Integrator.hpp
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/Integrator.h"
#include <cmath>

namespace gslpp {
namespace integration {

template<typename TArg, typename TRes, class Function>
void Integrator<TArg,TRes,Function>::integrate(TArg lborder, TArg uborder, Function &f,
		TRes &integral, gslpp::auxillary::NumAccuracyControl &integralAcc) {
	TRes errEstim;
	Gauss7Kronrad15(lborder,uborder,f,integral,errEstim);
}

template<typename TArg, typename TRes, class Function>
void Integrator<TArg,TRes,Function>::set_interval(TArg lborder,
		TArg uborder, Function const &f, IntervalG7K17 &intervalToBeSet) {

	TArg const kronradPoints[] = {-0.991455371120813,
			-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,
			-0.405845151377397,-0.207784955007898,0.000000000000000,0.207784955007898,
			0.405845151377397,0.586087235467691,0.741531185599394,0.864864423359769,
			0.949107912342759,0.991455371120813};

	//scale the points to the present interval
	for ( size_t i = 0 ; i < 15; ++i){
		intervalToBeSet.kronradPoints[i] = 0.5*( uborder*(kronradPoints[i]+1.0) - lborder*(kronradPoints[i]-1.0) );
	}

	for ( size_t i = 0 ; i < 15; ++i){

	}
	TRes const fEvals[] = {f(kronradPoints[0]),f(kronradPoints[1]),f(kronradPoints[2]),f(kronradPoints[3]),
			f(kronradPoints[4]),f(kronradPoints[5]),f(kronradPoints[6]),f(kronradPoints[7]),f(kronradPoints[8]),
			f(kronradPoints[9]),f(kronradPoints[10]),f(kronradPoints[11]),f(kronradPoints[12]),f(kronradPoints[13]),
			f(kronradPoints[14])};
}

template<typename TArg, typename TRes, class Function>
void  Integrator<TArg,TRes,Function>::s(TArg lborder, TArg uborder, Function &f,
		TRes &integral, TRes &errEstim) {

	//Gauss-Kronrad points and weights in the interval -1,1 according to Wikipedia
	//	http://en.wikipedia.org/wiki/Gauss%E2%80%93Kronrod_quadrature_formula

	TRes kronradWeights[] = {0.022935322010529,0.063092092629979,0.104790010322250,
			0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,
			0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,
			0.063092092629979,0.022935322010529};
	TRes gaussWeights[] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,
			0.381830050505119,0.279705391489277,0.129484966168870};

	TArg kronradPoints[] = {-0.991455371120813,
			-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,
			-0.405845151377397,-0.207784955007898,0.000000000000000,0.207784955007898,
			0.405845151377397,0.586087235467691,0.741531185599394,0.864864423359769,
			0.949107912342759,0.991455371120813};

	//scale the points to the present interval
	for ( size_t i = 0 ; i < 15; ++i){
		kronradPoints[i] = ( uborder*(kronradPoints[i]+1.0) - lborder*(kronradPoints[i]-1.0) ) / 2.0;
	}

	TArg const intervalLengthBy2 = (uborder - lborder) / 2.0;
	for ( size_t i = 0 ; i < 15; ++i)
		kronradWeights[i] *= intervalLengthBy2;
	for ( size_t i = 0 ; i < 7; ++i)
		gaussWeights[i] *= intervalLengthBy2;

	TRes const fEvals[] = {f(kronradPoints[0]),f(kronradPoints[1]),f(kronradPoints[2]),f(kronradPoints[3]),
			f(kronradPoints[4]),f(kronradPoints[5]),f(kronradPoints[6]),f(kronradPoints[7]),f(kronradPoints[8]),
			f(kronradPoints[9]),f(kronradPoints[10]),f(kronradPoints[11]),f(kronradPoints[12]),f(kronradPoints[13]),
			f(kronradPoints[14])};

	integral = 0;
	for ( size_t i = 0 ; i < 15; ++i)
		integral += fEvals[i]*kronradWeights[i];
	TRes integralGauss = 0;
	for ( size_t i = 0 ; i < 7; ++i)
		integralGauss += fEvals[2*i+1]*gaussWeights[i];

	errEstim = std::pow(200.0*std::fabs(integralGauss-integral),1.5);
};

} /* namespace integration */
} /* namespace gslpp */
