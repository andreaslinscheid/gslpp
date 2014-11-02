/*
 * Integrator.h
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_INTEGRATION_INTEGRATOR_H_
#define GSLPP_INTEGRATION_INTEGRATOR_H_

#include "gslpp/auxillary/NumAccuracyControl.h"
#include "gslpp/auxillary/FunctionTraits.h"
#include "gslpp/auxillary/has_function_signature.h"
#include <vector>
#include <type_traits>

namespace gslpp {
namespace integration {

/**
 * 	A class that provides integration tools for the definite integral.
 *
 *
 */
template<class Function, size_t indexT=0>
class Integrator {
public:

	/**	we deduce the result type from the Function */
	typedef typename gslpp::auxillary::FunctionTraits<Function>::result_type result_type;

	static_assert(indexT < gslpp::auxillary::FunctionTraits<Function>::nargs,
			"Second template parameter must be < the number of arguments to the operator()!");
	/**	we deduce the argument i type from the Function */
	typedef typename gslpp::auxillary::FunctionTraits<Function>::template arg<indexT>::type argument_type;

	/**
	 * Compute the adaptive integral using the template Function.
	 *
	 * @param lborder lower integral border.
	 * @param uborder upper integral border.
	 * @return The approximate integral \f$\int_{lborder}^{uborder} f(x)\rm{d}x\f$
	 */
	void integrate(argument_type lborder, argument_type uborder, Function const &f,
			result_type &integral,
			gslpp::auxillary::NumAccuracyControl<result_type> &integralAcc);

	/**
	 *
	 */
private:

	//The following construct identifies weight_type as the type specified by result_type::value_type
	//	or result_type if no such type exists.
	template<typename T = result_type>
	class result_type_trait {
		template<typename U>
		static auto check_type(typename U::value_type* dummy) -> typename U::value_type;

		template<typename U>
		static auto check_type(...) -> U;
	public:
		typedef decltype(check_type<T>(0)) value_type;
	};
	typedef typename result_type_trait<>::value_type weight_type;

	typedef struct {
		argument_type lborder;
		argument_type uborder;
		result_type integralVal;
		result_type errEstim;
		size_t subdiv;
	} Interval;

	void add_interval_points(argument_type lborder,argument_type uborder,std::vector<argument_type> &points) const;

	void evaluate_integral_formula_for_interval(size_t indexOfIntervalInData,
			weight_type intervalLength,
			std::vector<result_type> const& evaluatedPoints,
			result_type &integralOfInterval,
			result_type &errorEstiamteOfIntegral) const;

	weight_type distance_argument_types(Function const& f,argument_type v1, argument_type v2) const;

	void get_kronrad_points(argument_type lborder,argument_type uborder, argument_type (&kronradPoints)[15] ) const;

	void get_kronrad_weights(result_type (&kronradWeights)[15] ) const;
	void get_gauss_weights(result_type (&gaussWeights)[7] ) const;

	void evaluate_several_points(std::vector<argument_type> const &points,
			Function const &f,
			std::vector<result_type> &setOfEvaluatedPoints) const;
};

} /* namespace integration */
} /* namespace gslpp */

#include "gslpp/integration/src/Integrator.hpp"
#endif /* GSLPP_INTEGRATION_INTEGRATOR_H_ */
