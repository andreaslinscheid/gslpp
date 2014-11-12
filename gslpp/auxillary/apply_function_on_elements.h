/*
 * apply_function_on_elements.h
 *
 *  Created on: Nov 5, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_APPLY_FUNCTION_H_
#define GSLPP_AUXILLARY_APPLY_FUNCTION_H_

#include "gslpp/auxillary/FunctionTraits.h"
#include "gslpp/auxillary/has_iterator.h"

namespace gslpp{
namespace auxillary{

namespace delegate{
template<class method, bool resultHasMultipleElements>
struct apply_function_on_elements_impl;
}

/**
 * Apply the function or functor 'method' to all elements of obj with the arguments args.
 *
 * In the most simple case method needs to take an object of its return type as a first argument and
 * then the arguments args... that are passed.
 * Thus apply_function_on_elements<method>(obj,args...) compiles if obj = method(obj,args...) does.
 * We determine if method applies to many elements if obj has a type that defines an iterator.
 * Then we iterate through the elements and call *it=method(*it,args...) on each element.
 *
 * @param obj The object to which obj=method(obj,args...) or, if it defines an iterator, *it=method(*it,args...) is applied.
 * @param args The arguements passed to the method.
 */
template<class method,typename ...TArgs>
void apply_function_on_elements(typename FunctionTraits<method>::result_type obj,TArgs ...args) {
	delegate::apply_function_on_elements_impl<
		method,
		has_iterator<typename FunctionTraits<method>::result_type>::value
		>::call(obj,args...);
};

}
}

#include "gslpp/auxillary/src/apply_function_on_elements.hpp"
#endif /* GSLPP_AUXILLARY_APPLY_FUNCTION_H_ */
