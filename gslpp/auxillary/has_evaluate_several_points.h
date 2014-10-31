/*
 * has_evaluate_several_points.h
 *
 *  Created on: Oct 28, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_
#define GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_

#include <type_traits>

namespace gslpp {
namespace auxillary {

template<typename, typename T>
struct has_evaluate_several_points {
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

// specialization that does the checking

template<typename Obj, typename Ret, typename... Args>
struct has_evaluate_several_points<Obj, Ret(Args...)> {
private:
    template<typename T>
    static constexpr auto check(T*)
    	-> typename
    	std::is_same<
    			decltype(std::declval<T>().evaluate_several_points( std::declval<Args>()... )),
    			Ret
    			>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<Obj>(0)) type;

public:
    static constexpr bool value = type::value;
};


}
}
#endif /* GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_ */
