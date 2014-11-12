/*
 * FunctionTraits.hpp
 *
 *  Created on: Nov 5, 2014
 *      Author: alinsch
 */

#include "gslpp/auxillary/FunctionTraits.h"

namespace gslpp {
namespace auxillary {

//The follwoing simply fills the body of every function trait object. Based on the specified templates
//	the number of arguments and their type, as well as the return type is recorded.
#define STRUCT_BODY_FUNCTION_TRAITS 	                                          \
static const size_t nargs = sizeof...(Args);                                      \
                                                                                  \
typedef ResultType result_type;                                                   \
                                                                                  \
template <size_t i>                                                               \
struct arg{                                                                       \
    typedef typename std::tuple_element<i, std::tuple<Args...> >::type type;      \
}

template <class T, bool TIsClass>
struct FunctionTraits_impl { };

//The following specializations determine information about the operator() of a functor
template <class T>
struct FunctionTraits_impl<T,true> : public FunctionTraits_impl<decltype(&T::operator()),true> { };

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(Functor::*)(Args...) const,true> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType const&(Functor::*)(Args...) const,true> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (Functor::*)(Args...),true> {
	STRUCT_BODY_FUNCTION_TRAITS;
};


//The following specialization determine information about function pointers
template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType* (*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

//The following specialization determine information about member function pointers
template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType(T::*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType& (T::*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

template<class T,typename ResultType, typename ...Args>
struct FunctionTraits_impl<ResultType* (T::*)(Args...),false> {
	STRUCT_BODY_FUNCTION_TRAITS;
};

#undef STRUCT_BODY_FUNCTION_TRAITS

} /* namespace auxillary */
} /* namespace gslpp */
