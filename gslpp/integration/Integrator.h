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
#include <vector>

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
	typedef typename gslpp::auxillary::FunctionTraits<Function>::result_type TRes;

	/**	we deduce the argument i type from the Function */
	typedef typename gslpp::auxillary::FunctionTraits<Function>::template arg<indexT>::type TArg;

	/**
	 * Compute the adaptive integral using the template Function.
	 *
	 * @param lborder lower integral border.
	 * @param uborder upper integral border.
	 * @return The approximate integral \f$\int_{lborder}^{uborder} f(x)\rm{d}x\f$
	 */
	void integrate(TArg lborder, TArg uborder, Function const &f,
			TRes &integral, gslpp::auxillary::NumAccuracyControl<TRes> &integralAcc);

	/**
	 *
	 */
private:

	typedef struct {
		TArg lborder;
		TArg uborder;
		TRes integralVal;
		TRes errEstim;
	} Interval;

	typedef TArg Kronradpoints[15];
	typedef TArg Gausspoints[7];

	void get_kronrad_points(TArg lborder,TArg uborder, TArg (&kronradPoints)[15] ) const;

	void get_kronrad_weights(TRes (&kronradWeights)[15] ) const;
	void get_gauss_weights(TRes (&gaussWeights)[7] ) const;


	std::vector<TArg> _evaluationBuffer;
};

} /* namespace integration */
} /* namespace gslpp */

#include "gslpp/integration/src/Integrator.hpp"
#endif /* GSLPP_INTEGRATION_INTEGRATOR_H_ */
