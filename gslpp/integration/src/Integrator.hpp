/*
 * Integrator.hpp
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/Integrator.h"
#include <cmath>
#include <type_traits>

namespace gslpp {
namespace integration {

template<class Function,size_t indexT>
void Integrator<Function,indexT>::integrate(TArg lborder, TArg uborder, Function const &f,
		TRes &integral, auxillary::NumAccuracyControl<TRes> &integralAcc) {
	TRes errEstim;
	Integrator::Interval interval;
	interval.errEstim = 0;
	interval.integralVal = 0;
	interval.lborder = lborder;
	interval.uborder = uborder;
	interval.subdiv = 0;
	typename std::vector<Integrator::Interval> intervalsToBeDone( 1, interval );
	typename std::vector<Integrator::Interval> intervalsToBeDoneNextLoop;
	bool converged;
	std::vector<TArg> points;
	integral = 0;
	errEstim = 0;

	std::vector<Integrator::Interval> intervals;

	do {
		converged = true;
		size_t numIntervalsThisLoop = intervalsToBeDone.size();

		//collect all points that need to be evaluated
		typename std::vector<Integrator::Interval>::iterator it;
		for (  it = intervalsToBeDone.begin() ;it != intervalsToBeDone.end(); ++it){
			TArg newKronradPoints[15];
			this->get_kronrad_points(it->lborder,it->uborder,newKronradPoints);
			for ( size_t i = 0 ; i < 15; ++i)
				points.push_back(newKronradPoints[i]);
		}

		//evaluate all points
		std::vector<TRes> valAtPoints;
		valAtPoints.reserve(points.size());
		this->evaluate_several_points(points,f,valAtPoints);
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
			localContribution *= 0.5 * std::fabs( intervalsToBeDone[i].uborder - intervalsToBeDone[i].lborder);
			integralGauss *= 0.5 * std::fabs( intervalsToBeDone[i].uborder - intervalsToBeDone[i].lborder);

			TRes localErrEstim = std::pow(200.0*std::fabs(integralGauss-localContribution),1.5);

			bool thisIntervalConverged = integralAcc.locally_sufficient(localErrEstim,localContribution);
			converged = converged and thisIntervalConverged;

			if ( (not thisIntervalConverged) and integralAcc.sub_divisions_below_max(intervalsToBeDone[i].subdiv) ){
				TArg middle = 0.5 * ( intervalsToBeDone[i].uborder + intervalsToBeDone[i].lborder );
				intervalsToBeDone[i].subdiv += 1;
				Integrator::Interval intervalLower = intervalsToBeDone[i];
				Integrator::Interval intervalUpper = intervalsToBeDone[i];
				intervalLower.uborder = middle;
				intervalUpper.lborder = middle;
				intervalsToBeDoneNextLoop.push_back(intervalLower);
				intervalsToBeDoneNextLoop.push_back(intervalUpper);
			} else {
				errEstim += localErrEstim;
				integral += localContribution;
				intervals.push_back(intervalsToBeDone[i]);
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
				typename std::vector<Integrator::Interval>::iterator itMax = intervals.end();
				for ( it = intervals.begin(); it != intervals.end(); ++it){
					if ( maxErr < it->errEstim)
						itMax = it;
				}

				//the interval that should be subdivided is not dividable any more
				if ( integralAcc.sub_divisions_below_max(itMax->subdiv) )
					return;

				//remove the intervals contribution as the two subintervals will be re-added
				//in the next loop
				integral -= itMax->integralVal;
				errEstim -= itMax->errEstim;

				TArg middle = 0.5 * ( itMax->lborder + itMax->uborder );
				itMax->subdiv += 1;
				Integrator::Interval intervalLower = *itMax;
				Integrator::Interval intervalUpper = *itMax;
				intervalLower.uborder = middle;
				intervalUpper.lborder = middle;
				intervalsToBeDoneNextLoop.push_back(intervalLower);
				intervalsToBeDoneNextLoop.push_back(intervalUpper);
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

//We delegate to two implementations, one is simply evaluating the function using the operator()
//	and one is calling the function evaluate_several_points(points,setOfEvaluatedPoints) if the
//	object implements such a function (i.e. the Mode template parameter is true).
namespace delegate{

template<class Function,bool Mode>
struct evaluate_several_points_impl {};

//instantanated if Function _does_not_ implement the
//	evaluate_several_points(std::vector<TArg>,std::vector<TRes>) function
template<class Function>
struct evaluate_several_points_impl<Function,false> {
public:
	template<typename TRes,typename TArg>
	static void call (std::vector<TArg> const &points,
			Function const &f,
			std::vector<TRes> &setOfEvaluatedPoints){
		setOfEvaluatedPoints.clear();
		for (typename std::vector<TArg>::const_iterator it = points.begin(); it != points.end(); ++it){
			setOfEvaluatedPoints.push_back( f( *it) );
		}
	};
};

//instantanated if Function _does_ implement the
//	evaluate_several_points(std::vector<TArg>,std::vector<TRes>) function
template<class Function>
struct evaluate_several_points_impl<Function,true> {
public:
	template<typename TRes,typename TArg>
	static void call (std::vector<TArg> const &points,
			Function const &f,
			std::vector<TRes> &setOfEvaluatedPoints){
		f.evaluate_several_points(points,setOfEvaluatedPoints);
	};
};

};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::evaluate_several_points(std::vector<TArg> const &points,
		Function const &f,
		std::vector<TRes> &setOfEvaluatedPoints) const{

	//delegate to struct evaluate_several_points_impl above. See above documentation.
	delegate::evaluate_several_points_impl<Function,
		gslpp::auxillary::has_evaluate_several_points<Function, TRes(TArg) >::value
		>::template call<TRes,TArg>(points,f,setOfEvaluatedPoints);
};

} /* namespace integration */
} /* namespace gslpp */
