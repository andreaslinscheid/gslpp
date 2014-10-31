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

template <class Functor>
struct FunctionTraits : public FunctionTraits<decltype(&Functor::operator())>{
};

//specialize for pointer to const member function types
template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits<ResultType(Functor::*)(Args...) const>
{
    static const size_t nargs = sizeof...(Args);

    typedef ResultType result_type;

    template <size_t i>
    struct arg{
        typedef typename std::tuple_element<i, std::tuple<Args...> >::type type;
    };
};

//specialize for pointer to member function types
template<class Functor, typename ResultType, typename ...Args>
struct FunctionTraits<ResultType(Functor::*)(Args...)>
{
    static const size_t nargs = sizeof...(Args);

    typedef ResultType result_type;

    template <size_t i>
    struct arg{
        typedef typename std::tuple_element<i, std::tuple<Args...> >::type type;
    };
};

} /* namespace auxillary */
} /* namespace gslpp */

#endif /* GSLPP_AUXILLARY_FUNCTIONTRAITS_H_ */
