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

template<class Function,size_t indexT>
void Integrator<Function,indexT>::integrate(TArg lborder, TArg uborder, Function const &f,
		TRes &integral, auxillary::NumAccuracyControl<TRes> &integralAcc) {
	TRes errEstim;
	std::vector<std::pair<TArg,TArg> > intervalsToBeDone( 1, std::make_pair(lborder,uborder) );
	std::vector<std::pair<TArg,TArg> > intervalsToBeDoneNextLoop;
	bool converged;
	std::vector<TArg> points;
	integral = 0;
	errEstim = 0;

	std::vector<Integrator::Interval> intervals;

	do {
		converged = true;
		size_t numIntervalsThisLoop = intervalsToBeDone.size();

		//collect all points that need to be evaluated
		typename std::vector< std::pair<Integrator<Function,indexT>::TArg, Integrator<Function,indexT>::TArg> >::iterator it;
		for (  it = intervalsToBeDone.begin() ;it != intervalsToBeDone.end(); ++it){
			TArg newKronradPoints[15];
			this->get_kronrad_points(it->first,it->second,newKronradPoints);
			for ( size_t i = 0 ; i < 15; ++i)
				points.push_back(newKronradPoints[i]);
		}

		//evaluate all points
		std::vector<TRes> valAtPoints;
		for (typename std::vector<Integrator<Function,indexT>::TArg>::iterator it = points.begin(); it != points.end(); ++it){
			valAtPoints.push_back( f( *it) );
		}
		points.clear();

		//evaluate the integral and set up intervals for the next loop
		for ( size_t i = 0 ; i < numIntervalsThisLoop; i++) {

			TRes integralGauss = 0;
			TRes localContribution = 0;
			TRes kronradWeights[15];
			this->get_kronrad_weights(kronradWeights);
			TRes gaussWeights[7];
			this->get_gauss_weights(gaussWeights);
			for ( size_t point = 0 ; point < 15; ++point)
				localContribution += valAtPoints[i*15+point]*kronradWeights[point];
			for ( size_t point = 0 ; point < 7; ++point)
				integralGauss += valAtPoints[15*i+2*point+1]*gaussWeights[point];
			localContribution *= 0.5 * std::fabs( intervalsToBeDone[i].second - intervalsToBeDone[i].first);
			integralGauss *= 0.5 * std::fabs( intervalsToBeDone[i].second - intervalsToBeDone[i].first);

			TRes localErrEstim = std::pow(200.0*std::fabs(integralGauss-localContribution),1.5);

			bool thisIntervalConverged = integralAcc.locally_sufficient(localErrEstim,localContribution);
			converged = converged and thisIntervalConverged;

			if ( not thisIntervalConverged ){
				TArg middle = 0.5 * ( intervalsToBeDone[i].first + intervalsToBeDone[i].second );
				intervalsToBeDoneNextLoop.push_back(std::make_pair(intervalsToBeDone[i].first,middle));
				intervalsToBeDoneNextLoop.push_back(std::make_pair(middle,intervalsToBeDone[i].second));
			} else {
				errEstim += localErrEstim;
				integral += localContribution;
				Integrator::Interval interval;
				interval.errEstim = localErrEstim;
				interval.integralVal =localContribution;
				interval.lborder = intervalsToBeDone[i].first;
				interval.uborder = intervalsToBeDone[i].second;
				intervals.push_back(interval);
			}

		}

		//if local convergence has been achieved, check global convergence
		if ( converged ) {

			//at this point integral and errEstim represent the approximations for the
			//global integral and error estimation - check if this is sufficient
			if ( not integralAcc.global_sufficient(errEstim,integral)) {
				converged = false;

				//since global convergence has not been achieved,
				//subdivide the interval with the largest error estimate
				TRes maxErr = 0;
				typename std::vector<Integrator::Interval>::iterator it;
				typename std::vector<Integrator::Interval>::iterator itMax;
				for ( it = intervals.begin(); it != intervals.end(); ++it){
					if ( maxErr < it->errEstim )
						itMax = it;
				}

				//remove the intervals contribution as the two subintervals will be re-added
				//in the next loop
				integral -= itMax->integralVal;
				errEstim -= itMax->errEstim;

				TArg middle = 0.5 * ( itMax->lborder + itMax->uborder );
				intervalsToBeDoneNextLoop.push_back(std::make_pair(itMax->lborder ,middle));
				intervalsToBeDoneNextLoop.push_back(std::make_pair(middle,itMax->uborder));
			}
		}
		intervalsToBeDone.clear();
		intervalsToBeDone.swap(intervalsToBeDoneNextLoop);
	}while ( not converged );

}

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_kronrad_points(TArg lborder,TArg uborder, TArg (&kronradPoints)[15]) const{
	TArg tmp[15] = {-0.991455371120813,
			-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,
			-0.405845151377397,-0.207784955007898,0.000000000000000,0.207784955007898,
			0.405845151377397,0.586087235467691,0.741531185599394,0.864864423359769,
			0.949107912342759,0.991455371120813};
	//scale the points to the present interval

	for ( size_t i = 0 ; i < 15; ++i){
		kronradPoints[i] = 0.5*( uborder*(tmp[i]+1.0) - lborder*(tmp[i]-1.0) );
	};
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_kronrad_weights(TRes (&kronradWeights)[15] ) const {
	TRes tmp[15] = {0.022935322010529,0.063092092629979,0.104790010322250,
			0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,
			0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,
			0.063092092629979,0.022935322010529};

	for ( size_t i = 0 ; i < 15; ++i)
		kronradWeights[i] = tmp[i];
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_gauss_weights(TRes (&gaussWeights)[7] ) const {
	TRes tmp[7] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,
			0.381830050505119,0.279705391489277,0.129484966168870};

	for ( size_t i = 0 ; i < 7; ++i)
		gaussWeights[i] = tmp[i];
};

} /* namespace integration */
} /* namespace gslpp */
