/*
 * EvaluationMonitor.h
 *
 *  Created on: Oct 16, 2014
 *      Author: alinsch
 */

#ifndef GSLPP_AUXILLARY_EVALUATIONMONITOR_H_
#define GSLPP_AUXILLARY_EVALUATIONMONITOR_H_

#include <type_traits>

namespace gslpp {
namespace auxillary {

// Primary template with a static assertion
// for a meaningful error message
// if it ever gets instantiated.
// We could leave it undefined if we didn't care.

template<typename, typename T>
struct has_insert {
    static_assert(
        std::integral_constant<T, false>::value,
        "Second template parameter needs to be of function type.");
};

// specialization that does the checking

template<typename Obj, typename Ret, typename... Args>
struct has_insert<Obj, Ret(Args...)> {
private:
    template<typename T>
    static constexpr auto check(T*)
    	-> typename std::is_same<decltype( declval<T>().insert( std::declval<Args>()... ) ),Ret>::type;

    template<typename>
    static constexpr std::false_type check(...);

    typedef decltype(check<Obj>(0)) type;

public:
    static constexpr bool value = type::value;
};

template<class EvalMonitor, class Function>
struct monitor_evaluation {

	template<class U>
	static constexpr std::true_type test_has_store_function(
			decltype( &U::insert(Function::argument_type,Function::value_type) ) ) ;

	template<class U>
	static constexpr std::false_type test_has_store_function(...);

	static bool value = ( sizeof( test_has_store_function<EvalMonitor>(0) ) == sizeof(oneByte) );
};


template<typename TArg, typename TVal, class container>
class EvaluationMonitor {
};

} /* namespace auxillary */
} /* namespace gslpp */

#include "gslpp/auxillary/src/EvaluationMonitor.hpp"
#endif /* GSLPP_AUXILLARY_EVALUATIONMONITOR_H_ */
