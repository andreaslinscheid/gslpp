/*
 * has_function_signature.h
 *
 *  Created on: Nov 2, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_
#define GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_
/** \file has_function_signature.h
    \brief Checking if an object has a function of given signature.
*/

#include <type_traits>

namespace gslpp {
namespace auxillary {

/** Macro that generates the test for a given function name.
 *
 * 	generate a struct that perform the check if the first template paramter implements a
 * 	function of given name and (second template parameter) of given signature
 * 	using the SFINAE principle.
 */
#define MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME(name) 							                            \
																										\
template<typename, typename T>												                            \
struct has_##name {															                            \
	static_assert(															                            \
		std::integral_constant<T, false>::value,							                            \
		"Second template parameter needs to be of function type.");							            \
};																			                            \
                                                                                                        \
template<typename Obj, typename Ret, typename... Args>                                                  \
struct has_##name<Obj, Ret(Args...)> {                                            					    \
private:                                                                                                \
    template<typename T>                                                                                \
    static constexpr auto check(T*)                                                                     \
    	-> typename                                                                                     \
    	std::is_same<                                                                                   \
    			decltype(std::declval<T>().name( std::declval<Args>()... )),         					\
    			Ret                                                                                     \
    			>::type;                                                                                \
                                                                                                        \
    template<typename>                                                                                  \
    static constexpr std::false_type check(...);                                                        \
                                                                                                        \
    typedef decltype(check<Obj>(0)) type;                                                               \
                                                                                                        \
public:                                                                                                 \
                                                                                                        \
    static constexpr bool value = type::value;                                                          \
};

/** Create the check has_evaluate_several_points<T,type_return(ArgT1,...)> */
MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME(evaluate_several_points);

/** Create the check has_distance<T,type_return(ArgT1,...)> */
MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME(distance);

#undef MAKE_CHECK_HAS_FUNCTION_OF_THIS_NAME

}
}
#endif /* GSLPP_AUXILLARY_HAS_EVALUATE_SEVERAL_POINTS_H_ */
