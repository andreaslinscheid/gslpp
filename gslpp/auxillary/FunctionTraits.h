/*
 * FunctionTraits.h
 *
 *  Created on: Oct 21, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_FUNCTIONTRAITS_H_
#define GSLPP_AUXILLARY_FUNCTIONTRAITS_H_

#include <tuple>
#include <type_traits>

namespace gslpp {
namespace auxillary {

template <typename T, bool TIsClass>
struct FunctionTraits_impl;

/**
 * 	Determines the type traits of the function like type T which can also be a functor.
 *
 * 	Defines the result type as FunctionTraits<T>::result_type, the number of arguments
 * 	as FunctionTraits<T>::nargs and the type of the arguement number i as
 * 	FunctionTraits<T>::template arg<i>::type.
 */
template <typename T>
struct FunctionTraits : public FunctionTraits_impl<T,std::is_class<T>::value> { };

} /* namespace auxillary */
} /* namespace gslpp */

#include "gslpp/auxillary/src/FunctionTraits.hpp"
#endif /* GSLPP_AUXILLARY_FUNCTIONTRAITS_H_ */
