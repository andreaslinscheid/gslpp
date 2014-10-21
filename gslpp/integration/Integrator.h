/*
 * Integrator.h
 *
 *  Created on: Sep 23, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_INTEGRATION_INTEGRATOR_H_
#define GSLPP_INTEGRATION_INTEGRATOR_H_

#include "gslpp/auxillary/NumAccuracyControl.h"
#include <map>

namespace gslpp {
namespace integration {

/**
 * 	A class that provides integration tools for the definite integral.
 *
 *
 */
template<typename TArg, class TRes, class Function>
class Integrator {
public:

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
		TArg kronradPoints[17];
		TRes functionVals[17];
	} IntervalG7K17;

	TRes const _kronradWeights[] = {0.022935322010529,0.063092092629979,0.104790010322250,
			0.140653259715525,0.169004726639267,0.190350578064785,0.204432940075298,0.209482141084728,
			0.204432940075298,0.190350578064785,0.169004726639267,0.140653259715525,0.104790010322250,
			0.063092092629979,0.022935322010529};

	TRes const _gaussWeights[] = {0.129484966168870,0.279705391489277,0.381830050505119,0.417959183673469,
			0.381830050505119,0.279705391489277,0.129484966168870};


	std::vector<TArg> _evaluationBuffer;

	void set_interval(TArg lborder, TArg uborder, Function const &f, IntervalG7K17 &intervalToBeSet);

	void Gauss7Kronrad15(TArg lborder, TArg uborder, Function &f, TRes &integral, TRes &errEstim);

	void sub_div(IntervalG7K17 &intervalOrig,
			Function &f,
			IntervalG7K17 &intervalNewLower,
			IntervalG7K17 &intervalNewUpper );
};

} /* namespace integration */
} /* namespace gslpp */

#include "gslpp/integration/src/Integrator.hpp"
#endif /* GSLPP_INTEGRATION_INTEGRATOR_H_ */
