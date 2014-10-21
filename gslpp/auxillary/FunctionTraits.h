/*
 * FunctionTraits.h
 *
 *  Created on: Oct 21, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_FUNCTIONTRAITS_H_
#define GSLPP_AUXILLARY_FUNCTIONTRAITS_H_

#include <tuple>

namespace gslpp {
namespace auxillary {

}
template<typename T>
struct FunctionTraits{
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

template<typename ResultType, typename ...Args>
struct FunctionTraits<ResultType(Args...)>
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
