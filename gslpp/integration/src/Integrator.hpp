/*
 * Integrator.hpp
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/Integrator.h"
#include "gslpp/auxillary/has_function_signature.h"
#include "gslpp/auxillary/has_iterator.h"
#include <cmath>
#include <type_traits>
#include <complex>

namespace gslpp {
namespace integration {

template<class Function,size_t indexT>
void Integrator<Function,indexT>::integrate(
		argument_type lborder, argument_type uborder,
		Function const &f,
		result_type &integral,
		auxillary::NumAccuracyControl<result_type> &integralAcc) const {

	//set integral and error estimates to zero
	this->set_to_zero(integral);
	_zeroOfResultType = integral;
	result_type errEstim(_zeroOfResultType);

	bool converged;

	//set up the first initial interval
	Integrator::Interval interval;
	interval.errEstim = _zeroOfResultType;
	interval.integralVal = _zeroOfResultType;
	interval.lborder = lborder;
	interval.uborder = uborder;
	interval.subdiv = 0;
	typename std::vector<Integrator::Interval> intervalsToBeDone( 1, interval );
	typename std::vector<Integrator::Interval> intervalsToBeDoneNextLoop;

	std::vector<argument_type> points;

	//a collection of intervals that meet at least the local error bounds.
	std::vector<Integrator::Interval> intervals;

	do {
		converged = true;

		//collect all points that need to be evaluated
		typename std::vector<Integrator::Interval>::iterator it;
		for (  it = intervalsToBeDone.begin() ;it != intervalsToBeDone.end(); ++it)
			this->add_interval_points(it->lborder,it->uborder,points);

		//evaluate all points
		std::vector<result_type> valAtPoints;
		valAtPoints.reserve(points.size());
		this->evaluate_several_points(points,f,valAtPoints);
		points.clear();

		//evaluate the integral and set up intervals for the next loop
		for ( size_t i = 0 ; i < intervalsToBeDone.size(); i++) {

			result_type localContribution(_zeroOfResultType);
			result_type localErrEstim(_zeroOfResultType);
			weight_type intervalLength = this->distance_argument_types(
					f,intervalsToBeDone[i].uborder,intervalsToBeDone[i].lborder);
			this->evaluate_integral_formula_for_interval(i,intervalLength,valAtPoints,localContribution,localErrEstim);

			bool thisIntervalConverged = integralAcc.locally_sufficient(localErrEstim,localContribution);

			//split all intervals in two for the next loop that are not locally converged sufficiently.
			if ( (not thisIntervalConverged) and
				 (integralAcc.sub_divisions_below_max(intervalsToBeDone[i].subdiv)) ){
				argument_type middle =  ( intervalsToBeDone[i].uborder + intervalsToBeDone[i].lborder )
						* weight_type(0.5);
				intervalsToBeDone[i].subdiv += 1;
				Integrator::Interval intervalLower = intervalsToBeDone[i];
				Integrator::Interval intervalUpper = intervalsToBeDone[i];
				intervalLower.uborder = middle;
				intervalUpper.lborder = middle;
				intervalsToBeDoneNextLoop.push_back(intervalLower);
				intervalsToBeDoneNextLoop.push_back(intervalUpper);
			} else {
				errEstim = errEstim + localErrEstim;
				integral = integral + localContribution;
				intervals.push_back(intervalsToBeDone[i]);
			}

			//change converged to false in some interval is not converged
			converged = converged and thisIntervalConverged;
		}

		//if local convergence has been achieved, check global convergence
		if ( converged ) {

			//at this point integral and errEstim represent the approximations for the
			//global integral and error estimation - check if this is sufficient
			if ( not integralAcc.global_sufficient(errEstim,integral)) {
				converged = false;

				//since global convergence has not been achieved,
				//subdivide the interval with the largest error estimate
				typename std::vector<Integrator::Interval>::iterator it;
				typename std::vector<Integrator::Interval>::iterator itMax = intervals.begin();
				for ( it = intervals.begin(); it != intervals.end(); ++it){
					if ( integralAcc.first_lower_than_second(itMax->errEstim,it->errEstim) )
						itMax = it;
				}

				//the interval that should be subdivided is not dividable any more
				if ( integralAcc.sub_divisions_below_max(itMax->subdiv) ){
					//Warn
					return;
				}

				//remove the intervals contribution as the two subintervals will be re-added
				//in the next loop
				integral = integral + itMax->integralVal*static_cast<weight_type>(-1.0);
				errEstim = errEstim + itMax->errEstim*static_cast<weight_type>(-1.0);;

				//add the two intervals to be done in the next loop
				argument_type middle = ( itMax->lborder + itMax->uborder )*0.5;
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
void Integrator<Function,indexT>::non_adaptive_integral(
		std::vector<argument_type> segmentPoints,
		Function const &f,
		result_type &integral,
		result_type &errorEstimation) const {

	//set integral and error estimates to zero
	this->set_to_zero(integral);
	_zeroOfResultType = integral;
	errorEstimation = _zeroOfResultType;

	//numerical integral of zero measure
	if ( segmentPoints.size() < 2 )
		return;

	//collect all points that need to be evaluated
	std::vector<argument_type> points;
	for ( size_t i = 0 ; i < segmentPoints.size() - 1; ++i) {
		this->add_interval_points(segmentPoints[i],segmentPoints[i+1],points);
	}

	//evaluate all points
	std::vector<result_type> valueFAtPoints;
	valueFAtPoints.reserve(points.size());
	this->evaluate_several_points(points,f,valueFAtPoints);

	//compute the integral
	for ( size_t i = 0 ; i < segmentPoints.size() - 1; ++i) {

		weight_type intervalLength =
				this->distance_argument_types(f,segmentPoints[i],segmentPoints[i+1]);

		result_type localContribution(_zeroOfResultType);
		result_type localErrEstim(_zeroOfResultType);

		this->evaluate_integral_formula_for_interval(i,intervalLength,valueFAtPoints,localContribution,localErrEstim);

		integral = integral + localContribution;
		errorEstimation = errorEstimation + localErrEstim;
	}
}

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_kronrad_points(argument_type lborder,argument_type uborder, argument_type (&kronradPoints)[15]) const{
	weight_type const tmp[15] = {-0.991455371120813,
			-0.949107912342759,-0.864864423359769,-0.741531185599394,-0.586087235467691,
			-0.405845151377397,-0.207784955007898,0.000000000000000,0.207784955007898,
			0.405845151377397,0.586087235467691,0.741531185599394,0.864864423359769,
			0.949107912342759,0.991455371120813};
	//scale the points to the present interval
	for ( size_t i = 0 ; i < 15; ++i){
		kronradPoints[i] = ( uborder*(tmp[i]+1.0) - lborder*(tmp[i]-1.0) )*0.5;
	};
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_kronrad_weights(weight_type (&kronradWeights)[15] ) const {
	weight_type tmp[15] = {0.022935322010529,0.063092092629979,0.104790010322250,
			0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,
			0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,
			0.063092092629979,0.022935322010529};

	for ( size_t i = 0 ; i < 15; ++i)
		kronradWeights[i] = tmp[i];
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_gauss_weights(weight_type (&gaussWeights)[7] ) const {
	weight_type tmp[7] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,
			0.381830050505119,0.279705391489277,0.129484966168870};

	for ( size_t i = 0 ; i < 7; ++i)
		gaussWeights[i] = tmp[i];
};

//We delegate to two implementations, one is simply evaluating the function using the operator()
//	and one is calling the function evaluate_several_points(points,setOfEvaluatedPoints) if the
//	object implements such a function (i.e. the Mode template parameter is true).
namespace delegate{

template<class Function,typename result_type,typename argument_type,bool Mode>
struct evaluate_several_points_impl {};

//instantanated if Function _does_not_ implement the
//	evaluate_several_points(std::vector<argument_type>,std::vector<result_type>) function
template<class Function,typename result_type,typename argument_type>
struct evaluate_several_points_impl<Function,result_type,argument_type,false> {
	static void call (std::vector<argument_type> const &points,
			Function const &f,
			std::vector<result_type> &setOfEvaluatedPoints){
		setOfEvaluatedPoints.clear();
		for (typename std::vector<argument_type>::const_iterator it = points.begin(); it != points.end(); ++it){
			setOfEvaluatedPoints.push_back( f( *it) );
		}
	};
};

//instantanated if Function _does_ implement the
//	evaluate_several_points(std::vector<argument_type>,std::vector<result_type>) function
template<class Function,typename result_type,typename argument_type>
struct evaluate_several_points_impl<Function,result_type,argument_type,true> {
	static void call (std::vector<argument_type> const &points,
			Function const &f,
			std::vector<result_type> &setOfEvaluatedPoints){
		f.evaluate_several_points(points,setOfEvaluatedPoints);
	};
};

}; /* namespace delegate */

template<class Function,size_t indexT>
void Integrator<Function,indexT>::evaluate_several_points(std::vector<argument_type> const &points,
		Function const &f,
		std::vector<result_type> &setOfEvaluatedPoints) const{

	//delegate to struct evaluate_several_points_impl above. See above documentation.
	delegate::evaluate_several_points_impl<Function,result_type,argument_type,
		gslpp::auxillary::has_evaluate_several_points<Function, result_type(argument_type) >::value
		>::call(points,f,setOfEvaluatedPoints);
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::add_interval_points(argument_type lborder,
		argument_type uborder,std::vector<argument_type> &points) const {
	argument_type newKronradPoints[15];
	this->get_kronrad_points(lborder,uborder,newKronradPoints);
	for ( size_t i = 0 ; i < 15; ++i)
			points.push_back(newKronradPoints[i]);
}

template<class Function,size_t indexT>
void Integrator<Function,indexT>::evaluate_integral_formula_for_interval(
		size_t indexOfIntervalInData,
		weight_type intervalLength,
		std::vector<result_type> const& evaluatedPoints,
		result_type &integralOfInterval,
		result_type &errorEstiamteOfIntegral) const{

	// we assume integralOfInterval is set to zero
#ifdef DEBUG_BUILD
	//TODO Check that integralOfInterval is zero
#endif

	weight_type kronradWeights[15];
	this->get_kronrad_weights(kronradWeights);
	for ( size_t point = 0 ; point < 15; ++point)
		integralOfInterval = integralOfInterval +
			evaluatedPoints[indexOfIntervalInData*15+point]*kronradWeights[point];
	integralOfInterval = integralOfInterval * (intervalLength * static_cast<weight_type>(0.5));

	result_type integralGauss = integralOfInterval;
	this->set_to_zero(integralGauss);
	weight_type gaussWeights[7];
	this->get_gauss_weights(gaussWeights);
	for ( size_t point = 0 ; point < 7; ++point)
		integralGauss = integralGauss +
			evaluatedPoints[15*indexOfIntervalInData+2*point+1]*gaussWeights[point];
	integralGauss = integralGauss * (intervalLength * static_cast<weight_type>(0.5));

	errorEstiamteOfIntegral = this->evaluate_error( integralGauss, integralOfInterval);
}

//We delegate to two implementations, one is simply evaluating the function using the operator()
//	and one is calling the function evaluate_several_points(points,setOfEvaluatedPoints) if the
//	object implements such a function (i.e. the Mode template parameter is true).
namespace delegate{
template <class FuncT,class argumentT,class weightT, bool implmentsDistance>
struct distance_argument_types_impl{ };

template <class FuncT,class argumentT,class weightT>
struct distance_argument_types_impl<FuncT,argumentT,weightT,false> {
	static weightT call(
			FuncT const& f,argumentT v1,argumentT v2) {
		return std::fabs(v1 - v2);
	}
};

template <class FuncT,class argumentT,class weightT>
struct distance_argument_types_impl<FuncT,argumentT,weightT,true> {
	static weightT call(
			FuncT const& f,argumentT v1,argumentT v2) {
		return f.distance(v1,v2);
	}
};
}/* namespace delegate */

template<class Function,size_t indexT>
typename Integrator<Function,indexT>::weight_type
		Integrator<Function,indexT>::distance_argument_types(Function const& f,argument_type v1, argument_type v2) const{
	return delegate::distance_argument_types_impl<
				Function,argument_type,weight_type,
				auxillary::has_distance<Function,weight_type(argument_type,argument_type)>::value
			>::call(f,v1,v2);
}

template<class Function,size_t indexT>
typename Integrator<Function,indexT>::weight_type
Integrator<Function,indexT>::gauss_kronrad_err_est(weight_type estimGauss, weight_type estimKronrad) const {
	return std::pow(200.0*std::fabs(estimGauss-estimKronrad),1.5);
}

namespace delegate{
template <typename T,  bool THasIterator = auxillary::has_iterator<T>::value>
struct assign_zero { };

template <typename T>
struct assign_zero<T,false> {
	//any sensible floating point type will convert
	//this integer 0 to its exact zero representation
	static void call (T & toBeSetZero) {toBeSetZero = 0;};
};
template <typename T>
struct assign_zero<T,true> {
	static void call (T & toBeSetZero) {
		for ( auto &&element : toBeSetZero) {
			//any sensible floating point type will convert
			//this integer 0 to its exact zero representation
			element = 0;
		}
	};
};
}/* namespace delegate */

template<class Function,size_t indexT>
void Integrator<Function,indexT>::set_to_zero(result_type & toBeSetZero) const{
	delegate::assign_zero<result_type>::call(toBeSetZero);
};

namespace delegate{
template <class method,typename T,  bool THasIterator>
struct estimate_error_impl{ };

template <class method,typename T>
struct estimate_error_impl<method,T,false>{
	static T call(method const& f,T const&estimateMethod1, T const &estimateMethod2){
		return f(estimateMethod1,estimateMethod2);
	};
};

template <class method,typename T>
struct estimate_error_impl<method,std::complex<T>,false> {
	static std::complex<T> call(method const& f,std::complex<T> const&estimateMethod1, std::complex<T> const&estimateMethod2){
		return std::complex<T>(
				f(estimateMethod1.real(),estimateMethod2.real()),
				f(estimateMethod1.imag(),estimateMethod2.imag()));
	};
};

//specialization for std::complex is to perform the error estimation for real and imaginary parts independently
template <class method,typename T>
struct estimate_error_impl<method,T,true> {
	static T call(method const& f,T const& estim1, T const& estim2){
		T result = estim1;
		typename T::const_iterator it1,it2;
		typename T::iterator itr;
		it1 = estim1.begin();
		it2 = estim2.begin();
		itr = result.begin();
		for ( ; it1 != estim1.end(); ++it1){
			 ++it2;
			 ++itr;
			 *itr = f(*it1,*it2);
		}
		return result;
	};
};
};/* namespace delegate*/
template<class Function,size_t indexT>
typename Integrator<Function,indexT>::result_type
Integrator<Function,indexT>::evaluate_error(result_type const & estimate1, result_type const & estimate2) const{
	return delegate::estimate_error_impl<ErrorEstimationFunctor,result_type,
			auxillary::has_iterator<result_type>::value>::call(_errorEstimationFunctor,estimate1,estimate2);
}

template<class Function,size_t indexT>
template<typename T>
struct Integrator<Function,indexT>::result_type_trait<T,false>{
		typedef T value_type;
};

template<class Function,size_t indexT>
template<typename T>
struct Integrator<Function,indexT>::result_type_trait<T,true>{
	typedef typename T::value_type value_type;
};

} /* namespace integration */
} /* namespace gslpp */
