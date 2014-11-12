/*
 * Integrator.hpp
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#include "gslpp/integration/Integrator.h"
#include "gslpp/auxillary/apply_function_on_elements.h"
#include "gslpp/auxillary/has_function_signature.h"
#include <cmath>
#include <type_traits>
#include <complex>

namespace gslpp {
namespace integration {

template<class Function,size_t indexT>
void Integrator<Function,indexT>::integrate(argument_type lborder, argument_type uborder, Function const &f,
		result_type &integral, auxillary::NumAccuracyControl<result_type> &integralAcc) {

	//set integral and error estimates to zero
	result_type errEstim = this->set_to_zero();
	integral = this->set_to_zero();

	bool converged;

	//set up the first initial interval
	Integrator::Interval interval;
	interval.errEstim = this->set_to_zero();
	interval.integralVal = this->set_to_zero();
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

			result_type localContribution,localErrEstim;
			weight_type intervalLength = this->distance_argument_types(
					f,intervalsToBeDone[i].uborder,intervalsToBeDone[i].lborder);
			this->evaluate_integral_formula_for_interval(f,i,intervalLength,valAtPoints,localContribution,localErrEstim);

			bool thisIntervalConverged = integralAcc.locally_sufficient(localErrEstim,localContribution);

			//split all intervals in two for the next loop that are not locally converged sufficiently.
			if ( (not thisIntervalConverged) and
				 (integralAcc.sub_divisions_below_max(intervalsToBeDone[i].subdiv)) ){
				argument_type middle =  ( intervalsToBeDone[i].uborder + intervalsToBeDone[i].lborder ) * 0.5;
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
void Integrator<Function,indexT>::get_kronrad_weights(result_type (&kronradWeights)[15] ) const {
	weight_type tmp[15] = {0.022935322010529,0.063092092629979,0.104790010322250,
			0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,
			0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,
			0.063092092629979,0.022935322010529};

	for ( size_t i = 0 ; i < 15; ++i)
		kronradWeights[i] = tmp[i];
};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::get_gauss_weights(result_type (&gaussWeights)[7] ) const {
	weight_type tmp[7] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,
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
//	evaluate_several_points(std::vector<argument_type>,std::vector<result_type>) function
template<class Function>
struct evaluate_several_points_impl<Function,false> {
public:
	template<typename result_type,typename argument_type>
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
template<class Function>
struct evaluate_several_points_impl<Function,true> {
public:
	template<typename result_type,typename argument_type>
	static void call (std::vector<argument_type> const &points,
			Function const &f,
			std::vector<result_type> &setOfEvaluatedPoints){
		f.evaluate_several_points(points,setOfEvaluatedPoints);
	};
};

};

template<class Function,size_t indexT>
void Integrator<Function,indexT>::evaluate_several_points(std::vector<argument_type> const &points,
		Function const &f,
		std::vector<result_type> &setOfEvaluatedPoints) const{

	//delegate to struct evaluate_several_points_impl above. See above documentation.
	delegate::evaluate_several_points_impl<Function,
		gslpp::auxillary::has_evaluate_several_points<Function, result_type(argument_type) >::value
		>::template call<result_type,argument_type>(points,f,setOfEvaluatedPoints);
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
		Function const& f,
		size_t indexOfIntervalInData,
		weight_type intervalLength,
		std::vector<result_type> const& evaluatedPoints,
		result_type &integralOfInterval,
		result_type &errorEstiamteOfIntegral) const{

	integralOfInterval = this->set_to_zero();
	weight_type kronradWeights[15];
	this->get_kronrad_weights(kronradWeights);
	for ( size_t point = 0 ; point < 15; ++point)
		integralOfInterval = integralOfInterval +
			evaluatedPoints[indexOfIntervalInData*15+point]*kronradWeights[point];
	integralOfInterval = integralOfInterval * (intervalLength * static_cast<weight_type>(0.5));

	result_type integralGauss = this->set_to_zero();
	weight_type gaussWeights[7];
	this->get_gauss_weights(gaussWeights);
	for ( size_t point = 0 ; point < 7; ++point)
		integralGauss = integralGauss +
			evaluatedPoints[15*indexOfIntervalInData+2*point+1]*gaussWeights[point];
	integralGauss = integralGauss * (intervalLength * static_cast<weight_type>(0.5));

	errorEstiamteOfIntegral = this->set_to_zero();
	auxillary::apply_function_on_elements<
			decltype(&Integrator<Function,indexT>::gauss_kronrad_err_est)
				>(integralGauss,errorEstiamteOfIntegral,integralOfInterval);
}

//We delegate to two implementations, one is simply evaluating the function using the operator()
//	and one is calling the function evaluate_several_points(points,setOfEvaluatedPoints) if the
//	object implements such a function (i.e. the Mode template parameter is true).
namespace delegate{
template <typename Function,size_t indexT, typename weightT, bool implmentsDistance>
struct distance_argument_types_impl{ };

template <typename Function,size_t indexT, typename weightT>
struct distance_argument_types_impl<Function,indexT,weightT,false> {
	static weightT call(
			Function const& f,
			typename Integrator<Function,indexT>::argument_type v1,
			typename Integrator<Function,indexT>::argument_type v2) {
		return std::fabs(v1 - v2);
	}
};

template <typename Function,size_t indexT, typename weightT>
struct distance_argument_types_impl<Function,indexT,weightT,true> {
	static weightT call(
			Function const& f,
			typename Integrator<Function,indexT>::argument_type v1,
			typename Integrator<Function,indexT>::argument_type v2) {
		return f.distance(v1,v2);
	}
};
}

template<class Function,size_t indexT>
typename Integrator<Function,indexT>::weight_type
		Integrator<Function,indexT>::distance_argument_types(Function const& f,argument_type v1, argument_type v2) const{
	return delegate::distance_argument_types_impl<
			Function,
			indexT,
			weight_type,
			gslpp::auxillary::has_distance<Function,weight_type(argument_type,argument_type)>::value
			>::call(f,v1,v2);
}

namespace delegate{
template <typename result_type, typename weight_type>
struct estimate_error_impl{
	static result_type call(result_type &estimateMethod1, result_type &estimateMethod2){
		gslpp::auxillary::apply_function_on_elements<gauss_kronrad_error_estimation>(estimateMethod1+estimateMethod2*(-1.0));
		};
private:
	weight_type gauss_kronrad_error_estimation(weight_type diffGaussKronrad);
};

template <typename result_type>
struct estimate_error_impl<result_type,result_type> {
	static std::complex<result_type> call(std::complex<result_type> &estimateMethod1, std::complex<result_type> &estimateMethod2){
		return std::complex<result_type>(
				std::pow(200.0*std::fabs(estimateMethod1.real()-estimateMethod2.real()),1.5),
				std::pow(200.0*std::fabs(estimateMethod1.imag()-estimateMethod2.imag()),1.5));
	};
};

//specialization for std::complex is to perform the error estimation for real and imaginary parts independently
template <typename result_type>
struct estimate_error_impl<std::complex<result_type>,result_type> {
	static std::complex<result_type> call(std::complex<result_type> &estim1, std::complex<result_type> &estim2){
		return std::complex<result_type>(
				estimate_error_impl<result_type,result_type>::call(estim1.real(),estim2.real()),
				estimate_error_impl<result_type,result_type>::call(estim1.imag(),estim2.imag())
				);
	};
};
};/* namespace delegate*/

template<class Function,size_t indexT>
typename Integrator<Function,indexT>::weight_type
Integrator<Function,indexT>::gauss_kronrad_err_est(weight_type dummy, weight_type estimGauss, weight_type estimKronrad) const {
	return std::pow(200.0*std::fabs(estimGauss-estimKronrad),1.5);
}


template<class Function,size_t indexT>
typename Integrator<Function,indexT>::result_type
Integrator<Function,indexT>::set_to_zero() const{
	result_type zero = result_type()*0.0;
//	auto zeroFct = [] (weight_type dummy){
//		return static_cast<weight_type>(0.0);
//	};
//	auxillary::apply_function_on_elements< decltype( zeroFct ) >( zero, weight_type() );
	return zero;
};

} /* namespace integration */
} /* namespace gslpp */
